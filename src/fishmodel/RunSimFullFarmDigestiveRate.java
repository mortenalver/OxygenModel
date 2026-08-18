/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package fishmodel;

import fishmodel.enkf.AssimSettings;
import fishmodel.enkf.EnsembleKF;
import fishmodel.enkf.MpiHandler;
import fishmodel.enkf.Util;
import fishmodel.hydraulics.CurrentMagicFields;
import fishmodel.hydraulics.SimpleTankHydraulics;
import fishmodel.pellets.*;
import fishmodel.sim.InputDataNetcdf;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import save.SaveNetCDF;
import ucar.nc2.NetcdfFileWriteable;

import java.io.File;
import java.io.IOException;
import java.text.NumberFormat;
import java.text.SimpleDateFormat;
import java.util.*;

/**
 *
 * @author malv
 */
public class RunSimFullFarmDigestiveRate {

    public static final double HYPOXIA_THRESHOLD = 6;

    public static LinearInterpolator interpol = new LinearInterpolator();

    /**
     * Setup for digestive oxygen study
     */
    public static void main(String[] args) {

        NumberFormat nf1 = NumberFormat.getNumberInstance(Locale.ENGLISH);
        nf1.setMaximumFractionDigits(1);
        nf1.setMinimumFractionDigits(0);
        // Initialize random number generator:
        Random rnd = new Random();

        // Modelloppløsning:
        //double dxy = 2., dz = 2.; // Model resolution (m)
        double dxy = 3., dz = 3.; // Model resolution (m)

        boolean isControlSim = true;
        // Toggle whether we run the Bjørøya case or the artificial case
        boolean useBjoroyaData = true;

        boolean simulateStarving = false;

        boolean useNewConsumptionModel = true; // Use updated O2 model (PROHAV HI submitted 2025)

        // Save files:
        String saveDir = "./output_dig/";
        String simNamePrefix;
        if (useBjoroyaData) {
            // Bjørøya scenario:
            simNamePrefix = isControlSim ? "contrV6" : "digtestV6";
        } else {
            // Artificial scenario:
            //simNamePrefix = isControlSim ? "art_contr_2m" : "art_digtest_2m";
            simNamePrefix = isControlSim ? "tnew3_contr" : "tnew3_digtest";
        }
        if (simulateStarving)
            simNamePrefix = simNamePrefix + "_starved";


        String simNamePostfix = "";

        boolean doMPI = false; // Will be set to true if we are running is EnKF mode using MPI
        AssimSettings as = new AssimSettings(); // Settings related to the EnKF are gathered in AssimSettings
        if (as.dryRun)
            simNamePrefix += "_dr";

        int rank = 0, N=1;
        MpiHandler mpi = null;
        EnsembleKF enKF = null;
        try {
            mpi = new MpiHandler(args);
            System.out.println("rank="+mpi.getRank()+", N="+mpi.getN());
            rank = mpi.getRank();
            N = mpi.getN();
            doMPI = true;
            // TEST TEST TEST:
            AdvectPellets.disableMultiprocessing(); // No internal parallelization when running MPI
            //
            if (as.usePerturbations && ((rank < N-1) || !as.useTwin || as.perturbTwin))
                as.perturbThisMember = true;
        } catch (Throwable ex) {
            ex.printStackTrace();
            System.out.println("Not doing MPI");
        }
        boolean isRoot = (rank==0); // For convenience, isRoot tells us if this is the rank 0 process.

        boolean[][][] mask = null;
        double lastMeanFeedDepth = -1;
        boolean maskO2WhenSaving = false;

        boolean useConstantAmbientSaturation = false;
        double ambientSaturationVal = 1;

        boolean varyAmbient = false; // Reduction in ambient values towards the rest of the farm

        boolean useVerticalDist = true;

        boolean decreasingCurrentFactor = false;

        if (!useBjoroyaData)
            decreasingCurrentFactor = false; // Turn this off in artificial scenario

        // -----------------------------------------------------------
        // Activation of modified o2 uptake model:
        if (!isControlSim) {
            IngestionAndO2Subgrid.setAddDigestiveO2Cons(true); // If true, activating digestive o2 consumption
            if (useBjoroyaData) {

                if (!useNewConsumptionModel)
                    IngestionAndO2Subgrid.o2consumptionMultOld = 0.7882 * 1.3; // 3*1.3;
                else
                    IngestionAndO2Subgrid.o2consumptionMultNew = 0.65 * 1.0; // ??????
            }

            else {
                if (!useNewConsumptionModel)
                    IngestionAndO2Subgrid.o2consumptionMultOld = 0.95 * 0.7882 * 1.3; // Larger fish, higher V, need to down-adjust more to get equal means
                else
                    IngestionAndO2Subgrid.o2consumptionMultNew = 0.6791;
            }
        }
        // -----------------------------------------------------------

        int ndays = 1;

        // -----------------------------------------------------------

        double artifExtO2 = 8;
        if (!useBjoroyaData) {
            ndays = 16;
        }
        // -----------------------------------------------------------

        boolean useCurrentMagic = false; // Use spatially variable current flow field
        CurrentMagicFields cmf = null;
        if (useCurrentMagic) {
            cmf = new CurrentMagicFields("C:/Users/alver/OneDrive - NTNU/prosjekt/PROHAV/matlab/currents_heuristic_10deg.nc");
        }

        boolean includeHypoxiaAvoidance = true;
        int checkAvoidanceInterval = 30, checkAvoidanceCount = 0;

        int daysToAdd = 0;
        if (args.length >= 1) {
            for (int i=0; i<args.length; i++) {
                if (args[i].startsWith("offset:")) {
                    String lastPart = args[i].substring("offset:".length()).trim();
                    try {
                        daysToAdd = Integer.parseInt(lastPart);
                        break;
                    } catch (NumberFormatException e) {
                        throw new RuntimeException(e);
                    }
                }
            }
        }

        int initYear = 2022, initMonth = Calendar.JUNE, initDate = 22, initHour = 0, initMin = 0, initSec = 0;
        initDate += daysToAdd;
        double t_end = ndays*24.*3600.;//1*24*3600; // Duration of simulation
        int nSim = 9; // Number of days to simulate (separate sims)
        int startAt = 0; // Set to >0 to skip one of more simulations, but count them in the sim numbering



        // Domain settings and farm layout:
        double frameSize = 90; // Rammefortøyning
        double outerPadding = 75; // Ekstra rom utenfor rammefortøyningene
        int[] feedStartEnd = null;

        // Location setup
        int[] cageGrid = null;
        int[][] cagePos = null;
        double farmRotation = 0;
        //double nFish, meanWeight, feedPerDay;
        double[] count, meanWeight, feedPerDay;
        double rad = 25;
        double depth = 25, totDepth = 25; // Cage size (m)

        int storeIntervalFeed = 7200, storeIntervalInfo = 60;

        if (useBjoroyaData) {
            cageGrid = new int[]{4, 2};
            cagePos = new int[][]{{0, 0}, {0, 1}, {1, 0}, {2, 0}, {2, 1}, {3, 0}}; // cage index 4 is "our" cage
            // Actual cage numbers are (in order): 1, 7, 2, 9, 8, 4
            farmRotation = 42; // Current directions should be rotated by -1 times this angle
            feedStartEnd = new int[]{27000, 63000};

            //nFish = 169821; // Estimated number of individuals in experimental period (source: FishTalk data)
            //meanWeight = 2869.5; // Estimated mean weight in experimental period (source: FishTalk data)

            double[] biomass = new double[]
                    {446703, 113668, 414781, 442289, 478087, 438519};
            count = new double[]
                    {166727, 161771, 175341, 164988, 169810, 163232};
            meanWeight = new double[cagePos.length];
            feedPerDay = new double[cagePos.length];
            for (int i = 0; i < count.length; i++) {
                meanWeight[i] = 1000*biomass[i]/count[i];
            }


            //feedPerDay = 2900.*1000; // Approximate feeding for the one cage over 10 hours based on FishTalk data
        }
        else { // Artificial scenario
            cageGrid = new int[]{3, 3};
            cagePos = new int[][]{{0, 0}, {0, 1}, {0, 2}, {1, 0}, {1, 1}, {1, 2},
                    {2, 0}, {2, 1}, {2, 2}}; // cage index 4 is "our" cage (here, too)
            feedStartEnd = new int[]{43200, 86400};

            count = new double[cagePos.length];
            meanWeight = new double[cagePos.length];
            feedPerDay = new double[cagePos.length];
            for (int i = 0; i < count.length; i++) {
                count[i] = 200000; // Max number per cage
                meanWeight[i] = 4500; // Gives biomass of a little less than 20 kg/m3

            }
            // REDUCED RESOLUTION FOR TESTING
            //dxy = 4;
            //dz = 4;
            // REDUCED TIME RESOLUTION FOR 3D FIELDS:
            storeIntervalFeed = 4*7200;
        }

        for (int i=0; i<count.length; i++) {
            System.out.println("Cage "+(i+1)+": count="+count[i]+", meanWeight="+meanWeight[i]);
        }

        // Sensor depths (all horizontal positions will be equipped with sensors at all depths:
        double[] sensorDepths = new double[] {5, 10, 15};
        // Angle positions of sensors at outer edge of each tank (0 degrees refers to north):
        double[] sensorAngles = new double[] {128.2948, 2.8445, 246.8427};

        double[] domainDims = new double[] {outerPadding + outerPadding + frameSize*cageGrid[0],
                2*outerPadding + frameSize*cageGrid[1]};
        System.out.println("Domain dims: "+domainDims[0]+" x "+domainDims[1]);
        ArrayList<double[]> cagePositions = new ArrayList<>();
        for (int i=0; i<cagePos.length; i++) {
            cagePositions.add(new double[] {outerPadding + frameSize*((double)(cagePos[i][0]) +0.5),
                outerPadding + frameSize*(((double)cagePos[i][1]) +0.5)});
            double[] pos = cagePositions.get(cagePositions.size()-1);
            System.out.println("Cage: "+pos[0]+" x "+pos[1]);

        }

        // Cage settings:
        double dt = .5 * dxy; // Time step (s)
        boolean storeO2Histograms = true;
        double depthDomain = 38;

        System.out.println("Resolution: "+dxy+" , "+dz);

        double fishMaxDepth = 38; // The maximum depth of the fish under non-feeding condition

        double currentReductionFactor = 0.45;//0.6;//0.52;//0.8; // Multiplier for inside current as function of outside
        // If we are using time-dependent current reduction factor, this value will be updated inside the nsim loop further down.

        // Environmental conditions:
        double currentSpeedInit = 2*0.04; // External current speed (m/s)
        double T_w = 14;


        // Oxygen diffusion constant. To be set dependent on current speed.
        double diffKappaO2 = 0.1, diffKappaO2Z = 0.1;

        int[] cageDims = new int[3];
        double[] currentOffset = new double[] {0,0,0};
        double[] currentOffset_r = new double[] {0,0,0}; // Perturbed global current vector


        // Wind speed (x, y components in m/s) affecting feed spreader:
        double[] windSpeed = new double[] {0, 0};

        // Pellet setup:
        double[] sizes = new double[] {3, 6, 9, 12};
        double[] speeds = new double[] {0.0773, 0.0815, 0.1284, 0.1421};
        int di = 2; // Index of chosen pellet size
        double pelletWeight=0.2; // Pellet weight (g)
        double sinkingSpeed = speeds[di];
        double kappa_ref = 0.00012;
        double kappa_add = 0.2;
        double refSize = 9;
        double diffKappa = kappa_ref*(kappa_add + Math.pow(sizes[di]/refSize,2));
        double kappa_z_mult = 25;
        double diffKappaZ = diffKappa*kappa_z_mult;


        // Set up cage dimensions and cage grid:
        cageDims[0] = (int)Math.ceil(domainDims[0]/dxy);
        cageDims[1] = (int)Math.ceil(domainDims[1]/dxy);
        cageDims[2] = (int)Math.ceil(depthDomain/dz)+1;
        mask = CageMasking.fullFarmMasking(cageDims, dxy, cagePositions, rad, false);
        boolean useWalls = false;

        System.out.println("Domain dimensions: ("+cageDims[0]+", "+cageDims[1]+", "+cageDims[2]+")");

        // Feeding setup:
        int[][] feedingPos = new int[cagePositions.size()][2];
        for (int i=0; i<cagePositions.size(); i++) {
            double[] cp = cagePositions.get(i);
            feedingPos[i][0] = (int)Math.round(cp[0]/dxy);
            feedingPos[i][1] = (int)Math.round(cp[1]/dxy);
            System.out.println("Feeding pos "+i+": "+feedingPos[i][0]+" / "+feedingPos[i][1]);
        }

        // Feeding periods (start/end in s):
        int[][] feedingPeriods = new int[ndays][2];
        for (int i=0; i<ndays; i++) {
            feedingPeriods[i][0] = 86400*i + feedStartEnd[0];
            feedingPeriods[i][1] = 86400*i + feedStartEnd[1];
        }
        int nPeriods = feedingPeriods.length;

        Object sourceTerm = null;

        double[] ambientValueFeed = new double[cageDims[2]];
        for (int i = 0; i < ambientValueFeed.length; i++) {
            ambientValueFeed[i] = 0;
        }

        // Set up array for ambient temperature:
        double[] ambientTemp = new double[cageDims[2]];
        for (int i = 0; i < ambientTemp.length; i++) {
            ambientTemp[i] = T_w;
        }

        // Oxygen sensor positions:
        Measurements.MeasurementSet ms = Measurements.setupSensorPositionsAllCages(dxy, dz, rad,
                sensorDepths, sensorAngles, farmRotation,
                false, frameSize, cagePositions);

        // Feed affinity:
        double[][][] affinity = new double[cageDims[0]][cageDims[1]][cageDims[2]];
        for (int i=0; i<cageDims[0]; i++)
            for (int j=0; j<cageDims[1]; j++)
                for (int k=0; k<cageDims[2]; k++) {
                    //double depth = (k+0.5)*dz;
                    //if (depth > 10) affinity[i][j][k] = Math.max(0, 1-(depth-10)/5);
                    //else if (depth < 2)
                    //   affinity[i][j][k] = 0.2;
                    //else
                    if (mask==null || mask[i][j][k])
                        affinity[i][j][k] = 1;

                }

        // O2 affinity:
        // Vertical distribution data based on telemetry (8 individuals):
        double[] affProfile_orig = new double[] {0.0110, 0.0913, 0.8601, 2.1406, 2.7774, 2.6903, 2.5195, 2.2987, 2.0137,
                1.7448, 1.5883, 1.3667, 1.2348, 1.0724, 0.9379, 0.7764, 0.7104, 0.5895, 0.5607, 0.4668, 0.3933,
                0.4009, 0.2935, 0.1801, 0.1260, 0.0787, 0.0457, 0.0304};

        double[] affProfile_flat = new double[] {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
        double[] affProfile = useVerticalDist ? affProfile_orig : affProfile_flat;


        double[] affDepths = new double[] {0.5000, 1.5000, 2.5000, 3.5000, 4.5000, 5.5000, 6.5000, 7.5000, 8.5000,
                9.5000, 10.5000, 11.5000, 12.5000, 13.5000, 14.5000, 15.5000, 16.5000, 17.5000, 18.5000, 19.5000,
                20.5000, 21.5000, 22.5000, 23.5000, 24.5000, 25.5000, 26.5000, 27.5000};
        /*double[] affProfile = new double[] {1, 1};
        double[] affDepths = new double[] {0, 30};*/
        double[] affinityProfile = new double[cageDims[2]];
        interpolateVertical(affinityProfile, affDepths, affProfile, cageDims[2], dz);

        /*for (int i = 0; i < affinityProfile.length; i++) {
            double v = affinityProfile[i];
            System.out.println("Affinity: "+v);
        }*/
        double[][][] o2Affinity = new double[cageDims[0]][cageDims[1]][cageDims[2]];
        double o2AffSum = setO2AffinityWithVerticalProfile(cageDims, dz, fishMaxDepth, mask, affinityProfile, o2Affinity, affinity);
        double[] o2AffSums = MultiCageUtils.getO2AffSums(cagePositions, dxy, rad, o2Affinity);
        //double o2AffSum = setO2AffinityWithVerticalProfileAndEdgeDecrease(cageDims, rad, dxy, dz, fishMaxDepth, mask, affinityProfile, o2Affinity, affinity);
        int availableCellsForO2Uptake = countAvailableCellsForOxygenUptake(cageDims, dz, fishMaxDepth, mask);

        // Oxygen
        double[] ambientValueO2 = new double[cageDims[2]];
        double[] ambientValueO2_r = new double[cageDims[2]]; // Possibly perturbed ambient value of O2
        for (int i = 0; i < ambientValueO2.length; i++) {
            ambientValueO2[i] = artifExtO2;
        }

        // If we are running in EnKF mode, let rank 0 initialize the EnKF class:
        if (doMPI && (rank==0)) {
            enKF = new EnsembleKF(simNamePrefix, cageDims, as.nPar, dxy, ms);
        }

        // Set up initial perturbation and parameter values:
        double ambientO2_perturb = 0;
        double[] current_perturb = new double[2];
        double o2Cons_perturb = 0; // Relative perturbation to total oxygen consumption to be updated per time step
        double o2Cons_perturb_r = 0; // The perturbation to apply at this particular time step - may be set
        // equal to o2Cons_perturb, or to the sum of o2Cons_perturb and an estimated consumption parameter.
        double[] parVal = new double[as.nPar];
   
        // --------------------------------------------------------------------------
        // Simulations to run
        // --------------------------------------------------------------------------

        for (int sim=0; sim<nSim; sim++) {
            if (sim<startAt)
                continue;
            int feedingPeriodPiv = 0;
            boolean isFeeding = false;

            Calendar c = Calendar.getInstance();
            c.set(initYear, initMonth, initDate, initHour, initMin, initSec);
            c.add(Calendar.DATE, sim);
            Date startTime = c.getTime();
            int dayDelta = c.get(Calendar.DAY_OF_MONTH)-22; //sim;//initDate - 22.;
            System.out.println("Day delta = "+dayDelta);

            // Set up correct feeding values for this day:
            if (useBjoroyaData) {

                if (decreasingCurrentFactor) {
                    currentReductionFactor = 0.8 - 0.1 - (dayDelta * 0.2/8.0);
                    //currentReductionFactor = 0.8 + 0.05 - (dayDelta * 0.2/8.0); // Original values
                    System.out.println("Current reduction factor: "+currentReductionFactor);
                }

                //  FishTalk estimated growth in the days June 22 to 30 for all cages:
                double[][] estBiomassGrowth = new double[][] {
                        {1709, 823, 0, 0, 3042, 2008, 0, 0, 0},
                        {753, 757, 783, 294, 717, 153, 823, 708, 782},
                        {1104, 456, 0, 0, 5339, 3835, 4150, 0, 0},
                        {1943, 837, 0, 0, 3198, 2228, 2540, 2607, 0},
                        {1173, 952, 839, 0, 2524, 2466, 1922, 2353, 2135},
                        {842, 776, 783, 0, 3033, 1925, 2238, 2477, 2460}};

                for (int i = 0; i < count.length; i++) {
                    feedPerDay[i] = estBiomassGrowth[i][dayDelta]*1000*1.18;
                }
                //feedPerDay = 2900.*1000; // Approximate feeding for the one cage over 10 hours based on FishTalk data

            }
            else { // Artificial scenario

                for (int i = 0; i < count.length; i++) {
                    feedPerDay[i] = count[i]*meanWeight[i]*0.006; // 0.6 % of biomass per day

                }

            }


            // Current field
            double[][][][] hydro;
            // Here you set up the current profile (3D current vector per depth layer):
            double[][] currentProfile = new double[cageDims[2] + 1][3];
            for (int i=0; i<cageDims[2]+1; i++) {
                currentProfile[i][0] = 0;
            }
            hydro = SimpleTankHydraulics.getProfileHydraulicField(cageDims, currentProfile);


            AdvectPellets ap = new AdvectPellets();

            AdvectPellets apOx = new AdvectPellets();
            //AdvectPelletsVarCurr apOx_vc = new AdvectPelletsVarCurr();
            if (varyAmbient) {
                apOx.setVaryAmbient(true, affinityProfile);
                //apOx_vc.setVaryAmbient(true, affinityProfile);
            }

            // Format a unit string for the time variable to save to NetCDF giving the initial time:
            SimpleDateFormat formatter = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
            String unitString = "seconds since "+formatter.format(startTime);
            System.out.println("Unit string: "+unitString);

            // Initialize environmental input data:
            InputDataNetcdf inData = null;
            if (useBjoroyaData) {
                String inDataFile = "C:/Users/alver/OneDrive - NTNU/prosjekt/O2_Bjørøya/bjoroya_data.nc";
                if (!(new File(inDataFile)).exists())
                    inDataFile = "bjoroya_data.nc";
                inData = new InputDataNetcdf(inDataFile, true);
                inData.setStartTime(startTime);
            }


            SimpleFish[] fish = new SimpleFish[cagePositions.size()];
            double totN = 0, totWeight = 0;
            for (int i=0; i<fish.length; i++) {

                // Initialize simple (grouped) fish model:
                fish[i] = new SimpleFish(count[i], meanWeight[i], 0.2*meanWeight[i]);
                totN += count[i];
                totWeight += count[i]*meanWeight[i];

                if (useBjoroyaData && !simulateStarving) {
                    // Set initial V values

                    // The following array contains the final average V values per cage for the 8 first simulation days,
                    // preceded by a copy of the values for the first simulation day. This allows initialization of V
                    // values corresponding approximately to what should be expected based on feeding on the previous day.
                    double[][] initValsV = {
                            {3.3108, 1.6243, 2.1188, 3.7710, 2.3032, 1.7758}, // At the moment just copy of end of 22nd june
                            {3.3108, 1.6243, 2.1188, 3.7710, 2.3032, 1.7758},
                            {1.6992, 1.6225, 0.9995, 1.7402, 1.9000, 1.6446},
                            {0.2140, 1.6379, 0.2140, 0.2140, 1.6673, 1.6249},
                            {0.1840, 0.6825, 0.1840, 0.1840, 0.1840, 0.1840},
                            {5.2985, 1.4294, 7.8667, 5.6173, 4.3511, 5.3924},
                            {3.4557, 0.4220, 6.1426, 3.8550, 4.1334, 3.3872},
                            {0.1996, 1.6603, 6.9779, 4.6190, 3.4488, 4.1354},
                            {0.1535, 1.2327, 0.1535, 4.0499, 3.5705, 3.8954}
                    };
                    // Multiplier for each of the groups' V compared to the average:
                    double[] groupVMult = new double[] {0.4558, 0.6568, 0.8341, 1.0000, 1.1735, 1.3776, 1.6381};
                    // For Bjørøya comparison:
                    double[] initValCages = initValsV[sim];
                    double[] VtoSet = new double[7];
                    for (int j = 0; j < VtoSet.length; j++) {
                        VtoSet[j] = initValCages[i]*groupVMult[j];
                    }
                    fish[i].setAllV(VtoSet);
                            //new double[]{1.9702, 2.8394, 3.6058, 4.3228, 5.0730, 5.9551, 7.0810});
                }
            }
            System.out.println("Tot N: "+totN+" , Avg weight: "+(totWeight/totN));
            double[][][] fishTmp = new double[fish[0].getNGroups()][1][1];

            // Determine nominal feeding rate:
            double feedPeriodLength = feedingPeriods[0][1] - feedingPeriods[0][0];
            double[] nominalFeedingRate = new double[count.length];
            for (int i=0; i<count.length; i++) {
                nominalFeedingRate[i] = feedPerDay[i]/feedPeriodLength; // Approximate feeding over 10 hours based on FishTalk data
                if (simulateStarving)
                    nominalFeedingRate[i] = 1e-3;

            }


            double feedingRateMult = 0; // Set each timestep
            System.out.print("Feeding rate = ");
            for (int i=0; i<count.length; i++)
                System.out.print(nominalFeedingRate[i]+", ");
            System.out.println("");

            NumberFormat nf = NumberFormat.getNumberInstance(Locale.US);
            nf.setMaximumFractionDigits(2);
            nf.setMinimumFractionDigits(2);

            // Initialize states:
            double[][][] fc = new double[cageDims[0]][cageDims[1]][cageDims[2]];
            double[][][] o2 = new double[cageDims[0]][cageDims[1]][cageDims[2]];
            double[][][] ingDist = new double[cageDims[0]][cageDims[1]][cageDims[2]];
            double[][][] o2consDist = new double[cageDims[0]][cageDims[1]][cageDims[2]];
            //System.out.println("Initial ambient: "+inData.getO2Ambient5());
            //AdvectPellets.initField(o2, inData.getO2Ambient5());//avO2);
            AdvectPellets.initField(ingDist, 0);
            double outFlow = 0., outFlow_net = 0.;

            // Initialize O2 field based on first ambient values:
            // Initialize O2 field based on first ambient values:
            if (useBjoroyaData) {
                // Bjørøya scenario
                double[] ambVal = new double[]{inData.getO2Ambient5(), inData.getO2Ambient10(), inData.getO2Ambient15()};
                interpolateVertical(ambientValueO2, new double[]{5, 10, 15}, ambVal, cageDims[2], dz);
                for (int i = 0; i < cageDims[0]; i++)
                    for (int j = 0; j < cageDims[1]; j++)
                        for (int k = 0; k < cageDims[2]; k++) {
                            o2[i][j][k] = ambientValueO2[k];
                        }
            } else {
                // Artificial scenario
                for (int i = 0; i < cageDims[0]; i++)
                    for (int j = 0; j < cageDims[1]; j++)
                        for (int k = 0; k < cageDims[2]; k++) {
                            o2[i][j][k] = artifExtO2;
                        }
            }


            // Setup of surface feeding:
            double[][][] feedingRate = new double[cageDims[0]][cageDims[1]][cageDims[2]];
            double[][] fTemp = new double[cageDims[0]][cageDims[1]];

            double[][] surfFeed = new double[cageDims[0]][cageDims[1]];
            double totFeed = 0;
            for (int i=0; i<feedingPos.length; i++) {
                totFeed += nominalFeedingRate[i];
                PelletSpreaderModel.setPelletDist(fTemp, feedingPos[i][0], feedingPos[i][1], dxy, 0, 35, 0, true, windSpeed, null);
                // Add to total distribution:
                for (int ii=0; ii<surfFeed.length; ii++)
                    for (int jj = 0; jj < surfFeed[ii].length; jj++) {
                        surfFeed[ii][jj] = surfFeed[ii][jj] + nominalFeedingRate[i]*fTemp[ii][jj];
                    }
            }
            System.out.println("totFeed = "+totFeed);
            for (int i=0; i<surfFeed.length; i++)
                for (int j=0; j<surfFeed[i].length; j++) {
                    feedingRate[i][j][0] = surfFeed[i][j]/totFeed;//((double)feedingPos.length);
                }
            sourceTerm = feedingRate;


            SimpleDateFormat filenameForm = new SimpleDateFormat("yyyy_MM_dd");
            String filePrefix = simNamePrefix+"_"+filenameForm.format(startTime);

            // Establish file names to write data to:
            NetcdfFileWriteable ncfile = null;
            NetcdfFileWriteable fishfile = null;
            NetcdfFileWriteable histfile = null;
            String ncfilePath = saveDir + filePrefix + simNamePostfix + (doMPI ? "_"+String.format("%02d", rank) : "")+".nc";
            String fishfilePath = saveDir + filePrefix + simNamePostfix + (doMPI ? "_"+String.format("%02d", rank) : "")+"_fish.nc";
            String histFilePath = saveDir + filePrefix + simNamePostfix + (doMPI ? "_"+String.format("%02d", rank) : "")+"_hist.nc";
            boolean firstStore3d = true, firstStoreScalars = true, firstStoreHist = true;

            double totFeedAdded = 0;

            // Variables for current in artificial scenario:
            double currentDirection = 90;
            double currentSpeed = currentSpeedInit;
            double beta_currentDir = 0.001; // the value of beta is 1/tau where tau is the time to 37% correlation
            double sigma_currentDir = 15;//10;

            double t = 0;
            int n_steps = (int) (t_end / dt);
            long stime = System.currentTimeMillis();
            for (int i = 0; i < n_steps; i++) {

                //System.out.println("t = "+t);
                double tMin = t / 60;

                if (useBjoroyaData) {
                    if (inData.advance(t) || (i == 0)) {
                        //System.out.println("Updating environment: t = "+t);

                        double[] obsCurrentDepths = inData.getCurrentDepths();
                        double[] tempVal = new double[]{inData.getTemperature5(), inData.getTemperature10(), inData.getTemperature15()};
                        interpolateVertical(ambientTemp, new double[]{5, 10, 15}, tempVal, cageDims[2], dz);

                        double[] ambVal = new double[]{inData.getO2Ambient5(), inData.getO2Ambient10(), inData.getO2Ambient15()};
                        interpolateVertical(ambientValueO2, new double[]{5, 10, 15}, ambVal, cageDims[2], dz);

                        if (useConstantAmbientSaturation) {
                            // Set ambient O2 value based on a fraction of maximum saturation
                            // as a function of temperature. We use the surface temperature
                            // to choose the saturation level for all layers, because maximum
                            // saturation increases with increasing pressure, so it is only
                            // limiting at the surface.
                            for (int j = 0; j < ambientValueO2.length; j++) {
                                ambientValueO2[j] = ambientSaturationVal *
                                        OxygenSolubility.getOxygenSolubility(ambientTemp[0]);
                            }
                        }

                        double[] obsCurrentProfile = inData.getExtCurrentSpeedProfile();
                        double[] obsCurrentDirProfile = inData.getExtCurrentDirProfile();
                        double[] obsCurrentComp1 = new double[obsCurrentProfile.length],
                                obsCurrentComp2 = new double[obsCurrentProfile.length];

                        // TODO: Since the model domain is rotated we need to adjust the direction to compensate.
                        for (int j = 0; j < obsCurrentComp1.length; j++) {
                            double redFacHere = currentReductionFactor;
                            /* if (obsCurrentProfile[j] < 0.02)
                                redFacHere = 0.333*(1.0 + 2.0*currentReductionFactor);
                            else if (obsCurrentProfile[j] < 0.05)
                                redFacHere = 0.25*(1.0 + 3.0*currentReductionFactor);*/

                            obsCurrentComp1[j] = redFacHere *
                                    obsCurrentProfile[j] * Math.sin((obsCurrentDirProfile[j] - farmRotation) * Math.PI / 180.);
                            obsCurrentComp2[j] = redFacHere *
                                    obsCurrentProfile[j] * Math.cos((obsCurrentDirProfile[j] - farmRotation) * Math.PI / 180.);
                        }


                        double[] interpProfile1 = new double[cageDims[2]],
                                interpProfile2 = new double[cageDims[2]];
                        interpolateVertical(interpProfile1, obsCurrentDepths, obsCurrentComp1, cageDims[2], dz);
                        interpolateVertical(interpProfile2, obsCurrentDepths, obsCurrentComp2, cageDims[2], dz);


                        for (int j = 0; j < interpProfile1.length; j++) {
                            currentProfile[j][0] = interpProfile1[j];
                            currentProfile[j][1] = interpProfile2[j];
                            currentProfile[j][2] = 0.;
                        }

                        SimpleTankHydraulics.getProfileHydraulicField(hydro, cageDims, currentProfile);

                        currentOffset[0] = 0.;
                        currentOffset[1] = 0.;

                    }
                } else {
                    // Artificial scenario:
                    double day = Math.floor(t/86400.);
                    currentSpeed = 0.16-0.01*day;

                    currentDirection = Util.updateGaussMarkov(currentDirection, beta_currentDir, sigma_currentDir, dt, rnd);

                    currentOffset[0] = currentReductionFactor*currentSpeed*Math.cos(currentDirection*Math.PI/180.);
                    currentOffset[1] = currentReductionFactor*currentSpeed*Math.sin(currentDirection*Math.PI/180.);
                }

                diffKappaO2 = Math.min(0.5, 10*Math.pow(currentReductionFactor*0.06,2)); // Math.min(0.5, 10*Math.pow(currentReductionFactor*0.04,2));
                diffKappaO2Z = 5.0*0.1*diffKappaO2;
                //System.out.println("DiffKappa O2: "+diffKappaO2);

                // Update feeding rate depending on preset feeding periods:
                if (!isFeeding) { // Not already feeding. Check if we should start:
                    if ((feedingPeriodPiv < feedingPeriods.length) && (t >= feedingPeriods[feedingPeriodPiv][0])) {
                        isFeeding = true;
                        feedingRateMult = totFeed;
                    } else
                        feedingRateMult = 0;
                } else { // Already feeding. Check if we should stop:
                    if (t >= feedingPeriods[feedingPeriodPiv][1]) {
                        isFeeding = false;
                        feedingRateMult = 0;
                        feedingPeriodPiv++; // Update so we start looking for next period
                    } else
                        feedingRateMult = totFeed;
                }

                /*if (includeHypoxiaAvoidance) {
                    checkAvoidanceCount++;
                    if (checkAvoidanceCount == checkAvoidanceInterval) {
                        checkAvoidanceCount = 0;
                        o2AffSum = setO2AffinityWithAvoidance(cageDims, dz, fishMaxDepth, mask, o2, affinity, o2Affinity);
                        //System.out.println("Updating O2 affinity. sum="+o2AffSum);

                    }
                }*/


                totFeedAdded += dt * feedingRateMult;

                o2Cons_perturb_r = 0.;

                // Perturb if we are using MPI, except if we are using a twin, and this is the twin, and the
                // twin is not to be perturbed.
                if (doMPI && as.perturbThisMember) {

                    /*// Perturb anywhere: Repeat a given number of times:
                    for (int allstatesrep=0; allstatesrep<as.allStatesNRep*4; allstatesrep++) {
                        // Perturb all states randomly:
                        // Pick a random point and a perturbation, and let it drop off by r^2
                        int pt1 = (int) Math.floor(cageDims[0] * Math.random()),
                                pt2 = (int) Math.floor(cageDims[1] * Math.random()),
                                pt3 = (int) Math.floor(cageDims[2] * Math.random());
                        double perturbVal = as.allStatesStd * rnd.nextGaussian();
                        for (int ii = 0; ii < cageDims[0]; ii++)
                            for (int jj = 0; jj < cageDims[1]; jj++)
                                for (int kk = 0; kk < cageDims[2]; kk++) {
                                    double distance = (dxy * Math.sqrt((ii - pt1) * (ii - pt1) + (jj - pt2) * (jj - pt2) + (kk - pt3) * (kk - pt3)) *
                                            as.allStatesDistMultiplier) - as.allStatesMinDist;
                                    if (distance <= 1)
                                        o2[ii][jj][kk] += perturbVal;
                                    else
                                        o2[ii][jj][kk] += perturbVal / (distance * distance);
                                }
                    }*/

                    ambientO2_perturb = Util.updateGaussMarkov(ambientO2_perturb, as.ambientO2Beta, as.ambientO2Std, dt, rnd);
                    for (int j = 0; j < ambientValueO2.length; j++) {
                        ambientValueO2_r[j] = ambientValueO2[j] + ambientO2_perturb;
                    }
                    for (int j = 0; j < 2; j++) {
                        current_perturb[j] = Util.updateGaussMarkov(current_perturb[j], as.currentBeta, as.currentStd, dt, rnd);
                        currentOffset_r[j] = currentOffset[j] + current_perturb[j];
                    }
                    o2Cons_perturb = Util.getGaussValue(as.o2ConsStd, rnd);//Util.updateGaussMarkov(o2Cons_perturb, as.o2ConsBeta, as.o2ConsStd, dt, rnd);
                    o2Cons_perturb_r = o2Cons_perturb;
                } else {
                    // Copy ambientValueO2 and currentOffset without perturbations:
                    System.arraycopy(ambientValueO2, 0, ambientValueO2_r, 0, ambientValueO2.length);
                    System.arraycopy(currentOffset, 0, currentOffset_r, 0, currentOffset.length);
                }

                // Perturb parameters according to their std. setting if we have any:
                if (doMPI && as.perturbThisMember && as.nPar > 0) {

                    for (int j=0; j<as.nPar; j++) {
                        parVal[j] += as.parStd[j]*dt*rnd.nextGaussian();
                    }

                    /*
                    // The first three parameters are perturbations to ambient O2 at 5, 10 and 15 m. We need to
                    // calculate a linear interpolaton of these to all model depths before applying it:
                    double[] ambO2Par = new double[] {parVal[0], parVal[1], parVal[2]}; // Make array of amb O2 related parameters
                    double[] interpAmbO2Par = new double[cageDims[2]];
                    interpolateVertical(interpAmbO2Par, as.parDepths, ambO2Par, cageDims[2], dz);

                    // Apply parameter values to model:
                    // Param 0: offset to ambient O2 values:
                    for (int j=0; j<cageDims[2]; j++)
                        ambientValueO2_r[j] += interpAmbO2Par[j];*/
                    for (int j = 0; j < ambientValueO2.length; j++) {
                        ambientValueO2_r[j] += parVal[0];
                    }

                    /*// Parameter number 4 is additional perturbation to total o2 consumption:
                    o2Cons_perturb_r += parVal[3];*/
                }

                //System.out.println("Feedingratemult = "+feedingRateMult);
                double[] r = ap.step(dt, fc, dxy, dz, useWalls, mask, sinkingSpeed, diffKappa, diffKappaZ, 
                        hydro, currentOffset, sourceTerm, feedingRateMult, ambientValueFeed);
                outFlow = r[0]; // Feed lost from grid (not used)
                outFlow_net = r[1]; // Feed lost from the unmasked part of the grid (feed lost through side)


                double[] o2OutFlow = apOx.step(dt, o2, dxy, dz, useWalls, mask, 0, diffKappaO2, diffKappaO2Z,
                        hydro, currentOffset_r, feedingRate, 0, ambientValueO2_r);



                double totalIntake = 0, o2ConsumptionRate = 0;
                double rho = 0;
                for (int ii=0; ii<cagePositions.size(); ii++) {
                    int[][] ranges = MultiCageUtils.getRanges(cagePositions.get(ii), dxy, rad);

                    double[] res = IngestionAndO2Subgrid.calculateIngestion(dt, fc, o2, affinity, o2Affinity,
                            o2AffSums[ii], ranges, ingDist, o2consDist, dxy, dz, mask, pelletWeight, ambientTemp, fish[ii],
                            0, useNewConsumptionModel);
                    totalIntake += res[0];
                    rho += res[1];
                    o2ConsumptionRate += res[2];
                }
                rho /= (double)(cagePositions.size());

                t = t + dt;


                if (i>0 && ((t/((double)storeIntervalFeed) - Math.floor(t/(double)storeIntervalFeed)) < 1e-5)) {
                    double elapsed = (double) ((System.currentTimeMillis() - stime)) / 60000.;
                    double fractionCompleted = ((double) i) / ((double) n_steps);
                    double remaining = (elapsed / fractionCompleted) - elapsed;
                    System.out.println("t = " + nf1.format(t) + " - Estimated time to complete: " + nf1.format(remaining) + " minutes");

                    if (firstStore3d) {
                        firstStore3d = false;
                        ncfile = SaveNetCDF.initializeFile(ncfilePath, cageDims, dxy, dz, unitString, ms);
                        SaveNetCDF.createCageVariables(ncfile, "feed", "ingDist", "o2", "o2consDist");

                        // Make string describing cage layout:
                        StringBuilder sb2 = new StringBuilder();
                        for (double[] cp : cagePositions) {
                            sb2.append(cp[0]).append(",").append(cp[1]).append(";");
                        }
                        ncfile.addGlobalAttribute("cagePositions", sb2.toString());
                        ncfile.addGlobalAttribute("cageRad", rad);
                    }
                    else {
                        try {
                            ncfile = NetcdfFileWriteable.openExisting(ncfilePath);
                        } catch (IOException e) {
                            e.printStackTrace();
                        }
                    }

                    SaveNetCDF.saveCageVariable(ncfile, t, "feed", fc, mask, true);
                    SaveNetCDF.saveCageVariable(ncfile, t, "ingDist", ingDist, mask, false);
                    SaveNetCDF.saveCageVariable(ncfile, t, "o2", o2, (maskO2WhenSaving ? mask : null), false);
                    SaveNetCDF.saveCageVariable(ncfile, t, "o2consDist", o2consDist, mask, false);

                    try {
                        ncfile.close();
                    } catch (IOException e) {
                        e.printStackTrace();
                    }

                }

                if (i>0 && ((t/((double)storeIntervalInfo) - Math.floor(t/(double)storeIntervalInfo)) < 1e-5)) {
                //if ((t - Math.floor(t) < dt) && (Math.floor(t) % storeIntervalInfo == 0)) {

                    if (firstStoreScalars) {
                        firstStoreScalars = false;
                        fishfile = SaveNetCDF.initializeFile(fishfilePath, new int[]{fish[0].getNGroups(), 1, cageDims[2]}, 1, 1, unitString, ms);
                        SaveNetCDF.createProfileVariable(fishfile, "appetite", 0);
                        SaveNetCDF.createProfileVariable(fishfile, "ingested", 0);
                        SaveNetCDF.createScalarVariables(fishfile, "rho", "feedingRate", "o2ConsumptionRate",
                                "min_O2", "mean_O2",
                                "meanFeedDepth", "d_meanFeedDepth", "totIngRate", "totIngested", "totFeed",
                                "waste", "waste_net");
                        SaveNetCDF.createProfileVariable(fishfile, "ext_O2", 2); // dim=2 means along z dim
                        SaveNetCDF.createProfileVariable(fishfile, "temperature", 2);
                        SaveNetCDF.createProfileVariable(fishfile, "ext_currentU", 2);
                        SaveNetCDF.createProfileVariable(fishfile, "ext_currentV", 2);
                        SaveNetCDF.createScalarVariables(fishfile, ms.names);
                        for (int ii=0; ii<cagePositions.size(); ii++) {
                            SaveNetCDF.createScalarVariable(fishfile, "Cage_"+(ii+1)+"_V");
                            SaveNetCDF.createScalarVariable(fishfile, "Cage_"+(ii+1)+"_min");
                            SaveNetCDF.createScalarVariable(fishfile, "Cage_"+(ii+1)+"_perc5");
                            SaveNetCDF.createScalarVariable(fishfile, "Cage_"+(ii+1)+"_perc10");
                            SaveNetCDF.createScalarVariable(fishfile, "Cage_"+(ii+1)+"_mean");
                            SaveNetCDF.createScalarVariable(fishfile, "Cage_"+(ii+1)+"_fracHypoxia");
                        }
                    }
                    else {
                        try {
                            fishfile = NetcdfFileWriteable.openExisting(fishfilePath);
                        } catch (IOException e) {
                            e.printStackTrace();
                        }
                    }

                    if (storeO2Histograms) {
                        if (firstStoreHist) {
                            firstStoreHist = false;
                            histfile = SaveNetCDF.initializeHistogramFile(histFilePath, MultiCageUtils.getBinEdges(), cagePos.length, unitString);
                        }
                        else {
                            try {
                                histfile = NetcdfFileWriteable.openExisting(histFilePath);
                            } catch (IOException e) {
                                e.printStackTrace();
                            }
                        }
                    }

                    double[] groupArray = new double[fish[0].getNGroups()];

                    // Find cage with biomass to save for (where we only save one):
                    int cage = 0;
                    while ((cage<fish.length-1) && (fish[cage].getTotalN() == 0))
                        cage++;


                    SaveNetCDF.saveScalarVariable(fishfile, t, "waste", outFlow, true);
                    SaveNetCDF.saveScalarVariable(fishfile, t, "waste_net", outFlow_net, false);
                    
                    double totI = 0, totIngRate = 0;
                    for (int cageI=0; cageI<cagePositions.size(); cageI++) {
                        for (int j = 0; j < fishTmp.length; j++) {
                            groupArray[j] = fish[cageI].getIngested(j);
                            totI += fish[cageI].getN(j) * fishTmp[j][0][0];
                            totIngRate += fish[cageI].getN(j) * fish[0].getIngRate(j);
                        }
                    }
                    SaveNetCDF.saveProfileVariable(fishfile, t, "ingested", 0, groupArray, false);




                    for (int j = 0; j < groupArray.length; j++) {
                        groupArray[j] = fish[cage].getAppetite(j);
                    }
                    SaveNetCDF.saveProfileVariable(fishfile, t, "appetite", 0, groupArray, false);

                    SaveNetCDF.saveScalarVariable(fishfile, t, "totIngested", totI, false);

                    SaveNetCDF.saveScalarVariable(fishfile, t, "totIngRate", totIngRate, false);

                    double totalFeed = 0;
                    for (int ii = 0; ii < fc.length; ii++)
                        for (int j = 0; j < fc[0].length; j++)
                            for (int k = 0; k < fc[0][0].length; k++)
                                totalFeed += fc[ii][j][k];
                    SaveNetCDF.saveScalarVariable(fishfile, t, "totFeed", totalFeed, false);

                    double meanFeedDepth = AdvectPellets.getMeanFeedDepth(fc, dz, mask);
                    SaveNetCDF.saveScalarVariable(fishfile, t, "meanFeedDepth", meanFeedDepth, false);

                    double dMeanFeedDepth = (lastMeanFeedDepth > 0) ? (meanFeedDepth - lastMeanFeedDepth) / dt : 0;
                    SaveNetCDF.saveScalarVariable(fishfile, t, "d_meanFeedDepth", dMeanFeedDepth, false);

                    SaveNetCDF.saveScalarVariable(fishfile, t, "rho", rho, false);


                    SaveNetCDF.saveScalarVariable(fishfile, t, "o2ConsumptionRate", o2ConsumptionRate, false);

                    // Save temperature (input value):
                    SaveNetCDF.saveProfileVariable(fishfile, t, "temperature", 2, ambientTemp, false);

                    // Save external O2 (input value):
                    SaveNetCDF.saveProfileVariable(fishfile, t, "ext_O2", 2, ambientValueO2, false);


                    // Save external current speed and direction(input values):
                    double[] currentComp = new double[cageDims[2]];
                    for (int j = 0; j < currentComp.length; j++) {
                        currentComp[j] = currentProfile[j][0] + currentOffset_r[0];;
                    }
                    SaveNetCDF.saveProfileVariable(fishfile, t, "ext_currentU", 2, currentComp, false);
                    for (int j = 0; j < currentComp.length; j++) {
                        currentComp[j] = currentProfile[j][1] + currentOffset_r[1];;
                    }
                    SaveNetCDF.saveProfileVariable(fishfile, t, "ext_currentV", 2, currentComp, false);

                    // Save feeding rage (input value):
                    SaveNetCDF.saveScalarVariable(fishfile, t, "feedingRate", feedingRateMult, false);

                    // Save minimum O2 value:
                    //double[] values = minValueMeanAndFracHypoxia(o2, mask);
                    //SaveNetCDF.saveScalarVariable(fishfile, t, "min_O2", values[0], false);
                    //SaveNetCDF.saveScalarVariable(fishfile, t, "mean_O2", values[1], false);

                    //SaveNetCDF.saveScalarVariable(fishfile, t, "frac_hypoxia", values[2], false);

                    ArrayList<CageStats> cageStats = MultiCageUtils.getCageStats(o2, mask, cagePositions, rad, dxy,
                        HYPOXIA_THRESHOLD);
                    for (int ii=0; ii<cageStats.size(); ii++) {
                        CageStats st = cageStats.get(ii);
                        SaveNetCDF.saveScalarVariable(fishfile, t, "Cage_"+(ii+1)+"_min",
                                st.stats[0], false);
                        SaveNetCDF.saveScalarVariable(fishfile, t, "Cage_"+(ii+1)+"_perc5",
                                st.stats[1], false);
                        SaveNetCDF.saveScalarVariable(fishfile, t, "Cage_"+(ii+1)+"_perc10",
                                st.stats[2], false);
                        SaveNetCDF.saveScalarVariable(fishfile, t, "Cage_"+(ii+1)+"_mean",
                                st.stats[4], false);
                        SaveNetCDF.saveScalarVariable(fishfile, t, "Cage_"+(ii+1)+"_fracHypoxia",
                                st.stats[5], false);
                        // Save average gut content for cage:
                        SaveNetCDF.saveScalarVariable(fishfile, t, "Cage_"+(ii+1)+"_V",
                                fish[ii].getAverageV(), false);
                    }


                    // Save o2 level at sensor positions:
                    for (int j=0; j<ms.names.length; j++) {
                        // Store value at sensor position:
                        SaveNetCDF.saveScalarVariable(fishfile, t, ms.names[j], o2[ms.pos[j][0]][ms.pos[j][1]][ms.pos[j][2]], false);
                        /*// Get values in a neighbourhood of the sensor to calculate spatial variability:
                        double nearStd = getStdAround(o2, o2Pos[j]);
                        SaveNetCDF.saveScalarVariable(fishfile, t, o2Names[j]+"_std", nearStd, false);*/
                    }

                    lastMeanFeedDepth = meanFeedDepth;

                    try {
                        fishfile.close();
                    } catch (IOException e) {
                        e.printStackTrace();
                    }

                    if (storeO2Histograms) {

                        SaveNetCDF.saveHistograms(histfile, t, cageStats);

                        try {
                            histfile.close();
                        } catch (IOException e) {
                            e.printStackTrace();
                        }
                    }
                }


            }


            double totI = 0;
            for (int cageI=0; cageI<cagePositions.size(); cageI++) {
                for (int j = 0; j < fishTmp.length; j++) {
                    totI += fish[cageI].getN(j) * fish[cageI].getIngested(j);
                    if (Double.isNaN(fish[cageI].getN(j)))
                        System.out.println("fish.getN is NaN for j=" + j);
                    if (Double.isNaN(fish[cageI].getIngested(j)))
                        System.out.println("fish.getIngested is NaN for j=" + j);
                }
            }

            System.out.println("totFeedAdded = " + totFeedAdded);
            System.out.println("totI = " + totI);
            System.out.println("Feed wastage: " + nf.format(100 * (totFeedAdded - totI) / totFeedAdded) + " %");

            try {
                ncfile.close();
                fishfile.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }


    private static void setAllVal(double[][][][] h, double[] v) {
        for (int i=0; i<h.length; i++)
            for (int j=0; j<h[i].length; j++)
                for (int k=0; k<h[i][j].length; k++) {
                    for (int z=0; z<v.length; z++)
                        h[i][j][k][z] = v[z];
                }
    }

    public static double standardDev(double[] data) {
        // The mean average
        double mean = 0.0;
        for (int i=0; i<data.length; i++) {
            mean += data[i];
        }
        mean /= data.length;
        // The variance
        double variance = 0.0;
        for (int i=0; i<data.length; i++) {
            variance += Math.pow(data[i]-mean, 2);
        }
        variance /= (data.length-1);
        // Standard Deviation
        return Math.sqrt(variance);
    }


    private static double[] minValueMeanAndFracHypoxia(double[][][] h, boolean[][][] mask) {
        double minval = Double.MAX_VALUE;
        double meanVal = 0;
        int hypoCells = 0, nonHypoCells = 0, cageCells = 0;
        for (int i=0; i<h.length; i++)
            for (int j=0; j<h[i].length; j++)
                for (int k=0; k<h[i][j].length; k++) {
                    if (mask[i][j][k]) {
                        cageCells++;
                        meanVal += h[i][j][k];
                        if (h[i][j][k] < minval)
                            minval = h[i][j][k];
                        if (h[i][j][k] < HYPOXIA_THRESHOLD)
                            hypoCells++;
                        else
                            nonHypoCells++;
                    }
                }
        return new double[] {minval, meanVal/((double)cageCells),
                ((double)hypoCells)/((double)(nonHypoCells+hypoCells))};
    }

    private static int countAvailableCellsForOxygenUptake(int[] cageDims, double dz, double fishMaxDepth,
                                                          boolean[][][] mask) {
        int res = 0;
        for (int i=0; i<cageDims[0]; i++)
            for (int j=0; j<cageDims[1]; j++)
                for (int k=0; k<cageDims[2]; k++) {
                    double lDepth = ((double) k + 0.5) * dz;
                    if ((mask == null || mask[i][j][k]) && (lDepth < fishMaxDepth)) {
                        res++;
                    }
                }

        return res;
    }

    private static double setO2AffinityWithVerticalProfileAndEdgeDecrease(int[] cageDims, double rad, double dxy, double dz, double fishMaxDepth, boolean[][][] mask,
                                                           double[] affinityProfile, double[][][] affinity,
                                                           double[][][] o2Affinity) {
        double o2AffSum = 0;


        for (int i=0; i<cageDims[0]; i++)
            for (int j=0; j<cageDims[1]; j++) {

                double xDist = dxy*(double)(i - cageDims[0]/2);
                double yDist = dxy*(double)(j - cageDims[1]/2);
                double distFromCenter = Math.sqrt(xDist*xDist + yDist*yDist);

                for (int k = 0; k < cageDims[2]; k++) {
                    double lDepth = ((double) k + 0.5) * dz;
                    if ((mask == null || mask[i][j][k]) && (lDepth < fishMaxDepth)) {
                        // Check oxygen avoidance criterion:
                        affinity[i][j][k] = affinityProfile[k];
                        o2Affinity[i][j][k] = affinityProfile[k];

                    } else {
                        affinity[i][j][k] = 0;
                        o2Affinity[i][j][k] = 0;

                    }

                    // TEST TEST TEST
                    // Decrease affinity towards edge:
                    double edgeFactor = 1.;
                    if (distFromCenter/rad > 0.5) {
                        edgeFactor -= Math.pow((distFromCenter-0.5*rad)/(0.5*rad), 4);
                        affinity[i][j][k] *= edgeFactor;
                        o2Affinity[i][j][k] *= edgeFactor;
                    }

                    o2AffSum += o2Affinity[i][j][k];
                }
            }
        return o2AffSum;
    }

    private static double setO2AffinityWithVerticalProfile(int[] cageDims, double dz, double fishMaxDepth, boolean[][][] mask,
                                                           double[] affinityProfile, double[][][] affinity,
                                                           double[][][] o2Affinity) {
        double o2AffSum = 0;


        for (int i=0; i<cageDims[0]; i++)
            for (int j=0; j<cageDims[1]; j++)
                for (int k = 0; k < cageDims[2]; k++) {
                    double lDepth = ((double) k + 0.5) * dz;
                    if ((mask == null || mask[i][j][k]) && (lDepth < fishMaxDepth)) {
                        // Check oxygen avoidance criterion:
                        affinity[i][j][k] = affinityProfile[k];
                        o2Affinity[i][j][k] = affinityProfile[k];
                        //affinity[i][j][k] = 1.;//affinityProfile[k];
                        //o2Affinity[i][j][k] = 1.;//affinityProfile[k];

                    } else {
                        affinity[i][j][k] = 0;
                        o2Affinity[i][j][k] = 0;

                    }

                    o2AffSum += o2Affinity[i][j][k];
                }

        return o2AffSum;
    }

    private static double setO2AffinityWithAvoidance(int[] cageDims, double dz, double fishMaxDepth, boolean[][][] mask, double[][][] o2,
                                                     double[][][] affinity, double[][][] o2Affinity) {
        double o2AffSum = 0;

        double surfaceAvoidanceDepth = 4; // Range of (gradual) surface avoidance, m
        double avoidanceThresh = 4;

        for (int i=0; i<cageDims[0]; i++)
            for (int j=0; j<cageDims[1]; j++)
                for (int k=0; k<cageDims[2]; k++) {
                    double lDepth = ((double)k+0.5)*dz;
                    if ((mask==null || mask[i][j][k]) && (lDepth < fishMaxDepth)) {
                        // Check oxygen avoidance criterion:
                        if (o2[i][j][k] > avoidanceThresh) {
                            affinity[i][j][k] = 1;
                            o2Affinity[i][j][k] = 1;
                        }
                        else {
                            affinity[i][j][k] = o2[i][j][k]/avoidanceThresh; // Linear decrease below threshold
                            o2Affinity[i][j][k] = o2[i][j][k]/avoidanceThresh; // Linear decrease below threshold
                        }

                        // Check surface avoidance criterion:
                        /*if (lDepth < surfaceAvoidanceDepth) {
                            affinity[i][j][k] *= lDepth/surfaceAvoidanceDepth;
                            o2Affinity[i][j][k] *= lDepth/surfaceAvoidanceDepth;
                        }*/

                    }
                    o2AffSum += o2Affinity[i][j][k];
                }
        return o2AffSum;
    }


    private static double getStdAround(double[][][] field, int[] pos) {
        double[] values = new double[9];
        values[0] = field[pos[0]-1][pos[1]-1][pos[2]];
        values[1] = field[pos[0]-1][pos[1]][pos[2]];
        values[2] = field[pos[0]-1][pos[1]+1][pos[2]];
        values[3] = field[pos[0]][pos[1]-1][pos[2]];
        values[4] = field[pos[0]][pos[1]][pos[2]];
        values[5] = field[pos[0]][pos[1]+1][pos[2]];
        values[6] = field[pos[0]+1][pos[1]-1][pos[2]];
        values[7] = field[pos[0]+1][pos[1]][pos[2]];
        values[8] = field[pos[0]+1][pos[1]+1][pos[2]];

        return standardDev(values);
    }

    /**
     * Set up vertical profile for model grid based on values at set depths. Extrapolate beyond end values
     * Interpolates using linear interpolator from Apache Commons Math.
     * @param res The array to put interpolated values into
     * @param depths Depths at which values are given (increasing values)
     * @param values Values at given depths.
     * @param kmax Number of vertical layers
     * @param dz Vertical resolution
     * @return Interpolated/extrapolated profile
     */
    public static void interpolateVertical(double[] res, double[] depths, double[] values, int kmax, double dz) {

        double minDepth = depths[0], maxDepth = depths[depths.length-1],
                topValue = values[0], bottomValue = values[values.length-1];
        PolynomialSplineFunction interp = interpol.interpolate(depths, values);
        for (int i=0; i<res.length; i++) {
            double currDepth = ((double)i + 0.5)*dz;
            if (currDepth < minDepth)
                res[i] = topValue;
            else if (currDepth > maxDepth)
                res[i] = bottomValue;
            else
                res[i] = interp.value(currDepth);
        }


    }
}
