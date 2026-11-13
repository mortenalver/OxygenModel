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
import fishmodel.sim.InputDataLoenngrunnen;
import fishmodel.sim.InventoryNMBUStudy;
import save.SaveNetCDF;
import ucar.nc2.NetcdfFileWriteable;

import java.io.File;
import java.io.IOException;
import java.text.NumberFormat;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.Arrays;

/**
 *
 * @author malv
 */
public class RunSimCircularCages {


    /**
     * Setup for Lønngrunnen oxygen study i PROHAV
     *
     * Posisjon Lønngrunnen: 63°56.250′ 9°46.767′ (Fiskeridirektoratets kartverktøy
     * (63.9375 / 9.7795)
     */
    public static void main(String[] args) {

        NumberFormat nf1 = NumberFormat.getNumberInstance(Locale.ENGLISH);
        nf1.setMaximumFractionDigits(1);
        nf1.setMinimumFractionDigits(0);
        // Initialize random number generator:
        Random rnd = new Random();

        // ── Load configuration from JSON file ─────────────────────────────
        // Pass the config file as the first CLI argument (must end in .json).
        // Falls back to "loenngrunnen.json" in the working directory.
        String configPath = (args.length > 0 && args[0].endsWith(".json")) ? args[0] : "loenngrunnen.json";
        if (args.length > 0 && args[0].endsWith(".json"))
            args = Arrays.copyOfRange(args, 1, args.length);
        SimConfig cfg;
        try {
            cfg = new SimConfig(configPath);
            System.out.println("Loaded configuration from: " + configPath);
        } catch (IOException e) {
            throw new RuntimeException("Failed to load config file: " + configPath, e);
        }

        // Model resolution:
        double dxy = cfg.dxy, dz = cfg.dz;
        String resStr = String.valueOf(dxy);

        boolean useFieldInputs             = cfg.useFieldInputs;
        boolean useCurrentMagic            = cfg.useCurrentMagic;
        boolean useNewConsumptionModel     = cfg.useNewConsumptionModel;
        boolean useDigestiveModel          = cfg.useDigestiveModel;
        double  o2EvenFraction             = cfg.o2EvenFraction;
        double  extraConsumptionMultiplier = cfg.extraConsumptionMultiplier;
        double  currentReductionFactor     = useCurrentMagic ? 1.0 : cfg.currentReductionFactor; // Current magic overrides reduction factor and sets it to 1

        boolean testEastCurrAgainstFluent  = cfg.testEastCurrAgainstFluent;
        boolean cage6BiomassOnly           = cfg.cage6BiomassOnly;
        boolean dynamicCurrReductionFactor = cfg.dynamicCurrReductionFactor;

        int    nSim    = cfg.nSim;
        int    startAt = cfg.startAt;

        String outputPrefix     = cfg.outputPrefix;

        // Save files:
        String saveDir        = cfg.saveDir;
        String simNamePrefix  = outputPrefix + "_" + (useCurrentMagic ? "mag_" : "nomag_") + resStr + "m";
        String simNamePostfix = cfg.simNamePostfix;


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
            //ex.printStackTrace();
            System.out.println("Not doing MPI");
        }
        boolean isRoot = (rank==0); // For convenience, isRoot tells us if this is the rank 0 process.

        boolean[][][] mask = null;
        double lastMeanFeedDepth = -1;
        boolean maskO2WhenSaving = cfg.maskO2WhenSaving;
        boolean varyAmbient      = cfg.varyAmbient;
        boolean useVerticalDist  = cfg.useVerticalDist;


        CurrentMagicFields cmf = null;
        if (useCurrentMagic) {
            String cmPath = cfg.currentMagicDir + cfg.currentMagicFileName;
            if (!(new File(cmPath)).exists())
                cmPath = cfg.currentMagicFileName;
            cmf = new CurrentMagicFields(cmPath);
        }

        boolean includeHypoxiaAvoidance = cfg.includeHypoxiaAvoidance;
        int checkAvoidanceInterval = cfg.checkAvoidanceInterval, checkAvoidanceCount = 0;

        // Simulation start date from config (overrideable via CLI args year/month/day):
        int year = cfg.startYear, month = cfg.startMonth - 1, day = cfg.startDay;
        // month is zero-based in Calendar (January = 0)

        // If we are not running in MPI mode, we check the first argument whether it indicates a number of
        // days to add to the start date (if we are running in MPI mode the arguments actually contain
        // rank numbers and number of ranks, so can't be used to indicate days to add):
        if (!doMPI) {
            if (args.length >= 3) {
                try {
                    year = Integer.parseInt(args[0]);
                    month = Integer.parseInt(args[1])-1; // Month is zero-based, so subtract 1
                    day = Integer.parseInt(args[2]); // ... while day isn't

                    nSim = 1; // If start date is specified, simulate one day only
                } catch (NumberFormatException e) {
                    throw new RuntimeException(e);
                }
            }
        }

        // Simulation start time:
        int initYear = year, initMonth = month, initDate = day, initHour = 0, initMin = 0, initSec = 0;
        double t_end = 24*3600;//1*24*3600; // Duration of simulation

        // Domain settings and farm layout (from config):
        int[]    cageGrid     = new int[]  {cfg.cageGridRows, cfg.cageGridCols};
        int[][]  cagePos      = cfg.cagePositions;
        double   frameSize    = cfg.frameSize;
        double   outerPadding = cfg.outerPadding;
        double[] spaghettiRad   = cfg.spaghettiRad;
        double[] spaghettiDepth = cfg.spaghettiDepth;
        double   farmRotation   = cfg.farmRotation;
        int[]    feedStartEnd   = new int[] {cfg.feedStartSeconds, cfg.feedEndSeconds};

        // Sensor layout (from config):
        double[] sensorDepths = cfg.sensorDepths;
        double[] sensorAngles = cfg.sensorAngles;
        double   sensorrad    = cfg.sensorRad;

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
        double rad = spaghettiRad[0];
        double dt = .5 * dxy; // Time step (s)
        int     storeIntervalFeed = cfg.storeIntervalFeed;
        int     storeIntervalInfo = cfg.storeIntervalInfo;
        boolean storeO2Histograms = cfg.storeO2Histograms;
        double  depthDomain       = cfg.depthDomain;
        System.out.println("Resolution: "+dxy+" , "+dz);

        double fishMaxDepth     = cfg.fishMaxDepth;
        double currentSpeedInit = cfg.currentSpeedInit;
        double T_w              = cfg.T_w;
        double avO2             = cfg.avO2;


        // Oxygen diffusion constant. To be set dependent on current speed.
        double diffKappaO2, diffKappaO2Z;

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
        mask = CageMasking.fullFarmMaskingSpaghetti(cageDims, dxy, dz, cagePositions, spaghettiRad, spaghettiDepth, false);
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
        // Fra Eskil (Bjørøya): måltidene varte fra ca. kl. 07:30-17:30, i gjennomsnitt.
        int[][] feedingPeriods = new int[][] {{feedStartEnd[0], feedStartEnd[1]},
                {86400+feedStartEnd[0], 86400+feedStartEnd[1]},
                {2*86400+feedStartEnd[0], 2*86400+feedStartEnd[1]},
                {3*86400+feedStartEnd[0], 3*86400+feedStartEnd[1]},
                {4*86400+feedStartEnd[0], 4*86400+feedStartEnd[1]}};
        int nPeriods = feedingPeriods.length;
        /*int[][] feedingPeriods = new int[][] {{1*3600, 2*3600}, {3*3600, 4*3600}, {5*3600, 6*3600},
                {7*3600, 8*3600}, {9*3600, 10*3600}, {11*3600, 12*3600}, {13*3600, 14*3600}, {15*3600, 16*3600},
                {17*3600, 18*3600}, {19*3600, 20*3600}, {21*3600, 22*3600}, {23*3600, 24*3600}, {25*3600, 26*3600},
                {27*3600, 28*3600}, {29*3600, 30*3600}, {31*3600, 32*3600}};*/
        /*for (int i=0; i<nPeriods; i++) {
            System.out.println("Feeding period "+(i+1)+": "+feedingPeriods[i][0]+" to "+feedingPeriods[i][1]);
        }*/
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
        Measurements.MeasurementSet ms = Measurements.setupSensorPositionsLoenngrunnen(dxy, dz, sensorrad,
                sensorDepths, sensorAngles, farmRotation, cagePositions);

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

        // O2 affinity vertical profile (from config):
        double[] affProfile_orig = cfg.affinityValues;
        double[] affProfile_flat = new double[cfg.affinityValues.length];
        java.util.Arrays.fill(affProfile_flat, 1.0);
        double[] affProfile = useVerticalDist ? affProfile_orig : affProfile_flat;
        double[] affDepths  = cfg.affinityDepths;
        /*double[] affProfile = new double[] {1, 1};
        double[] affDepths = new double[] {0, 30};*/
        double[] affinityProfile = new double[cageDims[2]];
        SimUtils.interpolateVertical(affinityProfile, affDepths, affProfile, cageDims[2], dz,0);

        /*for (int i = 0; i < affinityProfile.length; i++) {
            double v = affinityProfile[i];
            System.out.println("Affinity: "+v);
        }*/
        double[][][] o2Affinity = new double[cageDims[0]][cageDims[1]][cageDims[2]];
        double o2AffSum = SimUtils.setO2AffinityWithVerticalProfile(cageDims, dz, fishMaxDepth, mask, affinityProfile, o2Affinity, affinity);
        double[] o2AffSums = MultiCageUtils.getO2AffSums(cagePositions, dxy, rad, o2Affinity);
        /*for (int i=0; i<o2AffSums.length; i++)
            System.out.println("Cage "+(i+1)+": affsum="+o2AffSums[i]);*/
        //double o2AffSum = SimUtils.setO2AffinityWithVerticalProfileAndEdgeDecrease(cageDims, rad, dxy, dz, fishMaxDepth, mask, affinityProfile, o2Affinity, affinity);
        int availableCellsForO2Uptake = SimUtils.countAvailableCellsForOxygenUptake(cageDims, dz, fishMaxDepth, mask);


        // Oxygen
        double[] ambientValueO2 = new double[cageDims[2]];
        double[] ambientValueO2_r = new double[cageDims[2]]; // Possibly perturbed ambient value of O2
        for (int i = 0; i < ambientValueO2.length; i++) {
            ambientValueO2[i] = avO2;
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


            // Current field
            double[][][][] hydro;
            // Here you set up the current profile (3D current vector per depth layer):
            double[][] currentProfile = new double[cageDims[2] + 1][3];
            for (int i=0; i<cageDims[2]+1; i++) {
                currentProfile[i][0] = 0;
            }
            hydro = SimpleTankHydraulics.getProfileHydraulicField(cageDims, currentProfile);


            AdvectPelletsVarCurr ap = new AdvectPelletsVarCurr();

            AdvectPelletsVarCurr apOx = new AdvectPelletsVarCurr();
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
            String inDataFile = cfg.inputDataDir + cfg.inDataFileName;
            InputDataLoenngrunnen inData = null;
            if (useFieldInputs) {
                if (!(new File(inDataFile)).exists())
                    inDataFile = cfg.inDataFileName;
                inData = new InputDataLoenngrunnen(inDataFile);
                inData.setStartTime(startTime);
            }

            // Read weight and count data for the simulation time from farmer's estimates:
            String invDataFile = cfg.inputDataDir + cfg.invDataFileName;
            if (!(new File(invDataFile)).exists())
                invDataFile = cfg.invDataFileName;
            InventoryNMBUStudy inv = new InventoryNMBUStudy(invDataFile, startTime);
            double[] weight = inv.getWeight(); // Average weight per cage
            double[] count = inv.getCount(); // Estimated count per cage
            double[] feed = inv.getFeed(); // Feed per cage (g per day)

            if (testEastCurrAgainstFluent) {
                double fishCount = 200000, avgWeight = 5000.0;
                for (int i = 0; i < feed.length; i++) {
                    count[i] = fishCount;
                    weight[i] = avgWeight;
                    feed[i] = 0.015*fishCount*avgWeight;
                }
            }

            SimpleFish[] fish = new SimpleFish[cagePositions.size()];
            double totN = 0, totWeight = 0;
            for (int i=0; i<fish.length; i++) {
                // Initialize simple (grouped) fish model:
                double cageCount = count[i];
                if (cage6BiomassOnly && (i != 5))
                    cageCount = 1e-6;
                fish[i] = new SimpleFish(count[i], weight[i], 0.2 * weight[i]);
                totN += cageCount;
                totWeight += count[i] * weight[i];

                System.out.println("Cage "+(i+1)+": count="+fish[i].getTotalN()+", weight="+fish[i].getW(3)+", feed="+feed[i]);
            }
            System.out.println("Tot N: "+totN+" , Avg weight: "+(totWeight/totN));
            double[][][] fishTmp = new double[fish[0].getNGroups()][1][1];

            // -----------------------------------------------------------
            // Set O2 even fraction for O2 consumption model:
            IngestionAndO2Subgrid.setO2EvenFraction(o2EvenFraction);
            // -----------------------------------------------------------

            // -----------------------------------------------------------
            // Activation of modified o2 uptake model:
            if (useDigestiveModel) {
                IngestionAndO2Subgrid.setAddDigestiveO2Cons(true); // If true, activating digestive o2 consumption
                if (!useNewConsumptionModel)
                    IngestionAndO2Subgrid.o2consumptionMultOld = 0.7882 * 1.3; // 3*1.3;
                else
                    IngestionAndO2Subgrid.o2consumptionMultNew *= 0.65 * 1.0; // ??????
                // Set initial V values:
                // Look up yesterday's inventory value:
                Calendar c2 = Calendar.getInstance();
                c2.set(initYear, initMonth, initDate, initHour, initMin, initSec);
                c2.add(Calendar.DATE, sim);
                c2.add(Calendar.DATE, -1);
                Date yesterday = c2.getTime();
                //System.out.println("YESTERDAY'S DATE = "+yesterday);
                InventoryNMBUStudy invTmp = new InventoryNMBUStudy(invDataFile, yesterday);
                double[] yFeed = invTmp.getFeed(), yCount = invTmp.getCount(), yWeight = invTmp.getWeight();
                double[] normVVal = new double[] {0.1407, 0.3791, 0.6705, 1.0000, 1.3850, 1.8483, 2.4247};
                for (int i=0; i<fish.length; i++) {
                    if (yCount[i] > 0) {
                        double relFeed = yFeed[i] / (yCount[i] * yWeight[i]);
                        // y = 1556*x - 0.845   (linear approximation of relationship between relative feeding rate and end value of V for 4th group)
                        double approxV = 1556*relFeed - 0.845;
                        //System.out.println("Cage "+(i+1)+" relFeed = "+relFeed+", approxV = "+approxV);
                        double[] vToSet = new double[7];
                        for (int j=0; j<7; j++)
                            vToSet[j] = normVVal[j]*approxV;
                        fish[i].setAllV(vToSet);
                    }
                }

                /*double[] initV = new double[] {1.2020, 3.2397, 5.7295, 8.5451, 11.8352, 15.7935, 20.7190};
                for (int i=0; i<fish.length; i++)
                    fish[i].setAllV(initV);*/
            }
            // -----------------------------------------------------------


            // -----------------------------------------------------------
            // Apply extra multiplier for o2 consumption:
            if (useNewConsumptionModel)
                IngestionAndO2Subgrid.o2consumptionMultNew = IngestionAndO2Subgrid.o2consumptionMultNew_orig*extraConsumptionMultiplier;
            else
                IngestionAndO2Subgrid.o2consumptionMultOld = IngestionAndO2Subgrid.o2consumptionMultOld_orig*extraConsumptionMultiplier;
            // -----------------------------------------------------------


            // Determine nominal feeding rate:
            double nominalFeedingRate = 0;//cagePositions.size()*2900.*1000/(10*3600); // Approximate feeding over 10 hours based on FishTalk data
            for (int ii=0; ii<feed.length; ii++) {
                nominalFeedingRate += feed[ii];
            }
            nominalFeedingRate /= (10*3600);

            double feedingRateMult = 0; // Set each timestep
            System.out.println("Feeding rate = "+nominalFeedingRate);

            NumberFormat nf = NumberFormat.getNumberInstance(Locale.US);
            nf.setMaximumFractionDigits(2);
            nf.setMinimumFractionDigits(2);

            // Initialize states:
            double[][][] fc = new double[cageDims[0]][cageDims[1]][cageDims[2]];
            double[][][] o2 = new double[cageDims[0]][cageDims[1]][cageDims[2]];
            double[][][] ingDist = new double[cageDims[0]][cageDims[1]][cageDims[2]];
            double[][][] o2consDist = new double[cageDims[0]][cageDims[1]][cageDims[2]];
            AdvectPellets.initField(o2, avO2);
            AdvectPellets.initField(ingDist, 0);
            double outFlow = 0., outFlow_net = 0.;

            // Initialize O2 field based on first ambient values:
            double[] ambVal;
            if (useFieldInputs) {
                ambVal = new double[]{inData.getO2Ambient5(), inData.getO2Ambient10(), inData.getO2Ambient15()};
            }
            else {
                ambVal = new double[]{avO2, avO2, avO2};
            }
            SimUtils.interpolateVertical(ambientValueO2, new double[] {5, 10, 15}, ambVal, cageDims[2], dz);
            for (int i=0; i<cageDims[0]; i++)
                for (int j=0; j<cageDims[1]; j++)
                    for (int k=0; k<cageDims[2]; k++) {
                        o2[i][j][k] = ambientValueO2[k];
                    }

            // Setup of surface feeding:
            double[][][] feedingRate = new double[cageDims[0]][cageDims[1]][cageDims[2]];
            double[][] fTemp = new double[cageDims[0]][cageDims[1]];

            double[][] surfFeed = new double[cageDims[0]][cageDims[1]];
            double totFeed = 0;
            for (int i=0; i<feedingPos.length; i++) {
                totFeed += feed[i];
                PelletSpreaderModel.setPelletDist(fTemp, feedingPos[i][0], feedingPos[i][1], dxy, 0, 35, 0, true, windSpeed, null);
                // Add to total distribution:
                for (int ii=0; ii<surfFeed.length; ii++)
                    for (int jj = 0; jj < surfFeed[ii].length; jj++) {
                        surfFeed[ii][jj] = surfFeed[ii][jj] + feed[i]*fTemp[ii][jj];
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

            double t = 0;
            int n_steps = (int) (t_end / dt);
            long stime = System.currentTimeMillis();
            for (int i = 0; i < n_steps; i++) {

                //System.out.println("t = "+t);
                double tMin = t / 60;

                if (useFieldInputs) {
                    if (inData.advance(t) || (i == 0)) {
                        //System.out.println("Updating environment: t = "+t);

                        double[] obsCurrentDepths = inData.getCurrentDepths();
                        /*for (int j = 0; j < obsCurrentDepths.length; j++) {
                            System.out.println("Obs depth layer "+j+": "+obsCurrentDepths[j]);
                        }*/

                        double[] tempVal = new double[]{ inData.getTemperature5(), inData.getTemperature10(), inData.getTemperature15()};
                        SimUtils.interpolateVertical(ambientTemp, new double[]{5, 10, 15}, tempVal, cageDims[2], dz);

                        ambVal = new double[]{inData.getO2Ambient5(), inData.getO2Ambient10(), inData.getO2Ambient15()};
                        SimUtils.interpolateVertical(ambientValueO2, new double[]{5, 10, 15}, ambVal, cageDims[2], dz);

                        double[] obsCurrentProfile = inData.getExtCurrentSpeedProfile();
                        double[] obsCurrentDirProfile = inData.getExtCurrentDirProfile();
                        double[] obsCurrentComp1 = new double[obsCurrentProfile.length],
                                obsCurrentComp2 = new double[obsCurrentProfile.length];

                        /*for (int j = 0; j < obsCurrentDirProfile.length; j++) {
                            System.out.println("Dir layer "+j+": "+obsCurrentDirProfile[j]+" speed="+obsCurrentProfile[j]);
                        }*/
                        //System.out.println("Temp 0: "+tempVal[0]+", Amb 0: "+ambVal[0]+", Speed 0: "+obsCurrentProfile[0]);


                        // TODO: Since the model domain is rotated we need to adjust the direction to compensate.
                        for (int j = 0; j < obsCurrentComp1.length; j++) {
                            double currentSpeedHere = currentReductionFactor*obsCurrentProfile[j];
                            if (dynamicCurrReductionFactor) {
                                currentSpeedHere = Math.max(currentReductionFactor*obsCurrentProfile[j],
                                        Math.min(obsCurrentProfile[j], 0.005 + 0.7*obsCurrentProfile[j]));
                            }
                            obsCurrentComp1[j] = currentSpeedHere * Math.sin((obsCurrentDirProfile[j] - farmRotation) * Math.PI / 180.);
                            obsCurrentComp2[j] = currentSpeedHere * Math.cos((obsCurrentDirProfile[j] - farmRotation) * Math.PI / 180.);
                        }
                        //System.out.println("u="+obsCurrentComp1[0]+" / v="+obsCurrentComp2[0]+" / origdir="+obsCurrentDirProfile[0]+", farmrotation="+farmRotation);
                        //for (int j = 0; j < obsCurrentComp1.length; j++) {
                        //    System.out.println("i="+j+", u="+obsCurrentComp1[j]+" / v="+obsCurrentComp2[j]);
                        //}

                        double[] interpProfile1 = new double[cageDims[2]],
                                interpProfile2 = new double[cageDims[2]];
                        SimUtils.interpolateVertical(interpProfile1, obsCurrentDepths, obsCurrentComp1, cageDims[2], dz);
                        SimUtils.interpolateVertical(interpProfile2, obsCurrentDepths, obsCurrentComp2, cageDims[2], dz);

                        /*for (int j = 0; j < interpProfile1.length; j++) {
                            System.out.println("Interpolated: i="+j+", u="+interpProfile1[j]+" / v="+interpProfile2[j]);
                        }*/

                        if (testEastCurrAgainstFluent) {
                            // Override environmental conditions for this test:
                            double testCurrSpeedU = 0.1, testCurrSpeedV = 0, ambO2val = 10, ambTempval = 14;
                            for (int j = 0; j < interpProfile1.length; j++) {
                                interpProfile1[j] = testCurrSpeedU;
                                interpProfile2[j] = testCurrSpeedV;
                                ambientValueO2[j] = ambO2val;
                                ambientTemp[j] = ambTempval;
                            }
                        }

                        for (int j = 0; j < interpProfile1.length; j++) {
                            //if (j==0)
                            //    System.out.println("Interpolated current speed "+(j)+": "+interpProfile1[j]+" / "+interpProfile2[j]);
                            currentProfile[j][0] = interpProfile1[j];
                            currentProfile[j][1] = interpProfile2[j];
                            currentProfile[j][2] = 0.;
                        }

                        if (!useCurrentMagic) {
                            // Update current field using the new profile:
                            SimpleTankHydraulics.getProfileHydraulicField(hydro, cageDims, currentProfile);
                        } else {
                            double[] lDirections = new double[cageDims[2]],
                                    lSpeeds = new double[cageDims[2]];
                            for (int j = 0; j < cageDims[2]; j++) {
                                lSpeeds[j] = Math.sqrt(interpProfile1[j] * interpProfile1[j] + interpProfile2[j] * interpProfile2[j]);
                                lDirections[j] = Math.atan2(interpProfile1[j], interpProfile2[j]) * 180. / Math.PI;

                            }
                            cmf.setCurrentField(hydro, lSpeeds, lDirections);
                        }
                    }
                }
                else { // Artificial setting, no field data
                    // Update variable current speed offset:
                    if (!useCurrentMagic) {
                        // Juse use the global offset:
                        currentOffset[0] = 0.1;//currentReductionFactor*currentSpeed*Math.cos(currentDirection*Math.PI/180.);
                        currentOffset[1] = 0.0;//currentReductionFactor*currentSpeed*Math.sin(currentDirection*Math.PI/180.);
                    } else {
                        currentOffset[0] = 0;
                        currentOffset[1] = 0;
                        //double ucomp = -0.5, vcomp = -0.035;
                        double angle = 90 + 2. * (Math.random() - 0.5);
                        double speed = 0.1;//Math.sqrt(ucomp*ucomp+vcomp*vcomp); //0.1;
                        double[] lDirections = new double[cageDims[2]],
                                lSpeeds = new double[cageDims[2]];
                        for (int j = 0; j < cageDims[2]; j++) {
                            lSpeeds[j] = speed;
                            lDirections[j] = angle;
                            //System.out.println("ANGLE: "+(90-Math.atan2(vcomp, ucomp)*180./Math.PI));
                        }
                        cmf.setCurrentField(hydro, lSpeeds, lDirections);
                    }
                }

                diffKappaO2 = 0.023;//Math.min(0.5, 10*Math.pow(currentReductionFactor*0.06,2)); // Math.min(0.5, 10*Math.pow(currentReductionFactor*0.04,2));
                diffKappaO2Z = diffKappaO2; // TEST TEST TEST TEST TEST //0.5*diffKappaO2;
                //diffKappaO2 *= 2; // TEST TEST TEST TEST TEST
                //System.out.println("DiffKappa O2: "+diffKappaO2);

                // Update feeding rate depending on preset feeding periods:
                if (!isFeeding) { // Not already feeding. Check if we should start:
                    if ((feedingPeriodPiv < feedingPeriods.length) && (t >= feedingPeriods[feedingPeriodPiv][0])) {
                        isFeeding = true;
                        feedingRateMult = nominalFeedingRate;
                    } else
                        feedingRateMult = 0;
                } else { // Already feeding. Check if we should stop:
                    if (t >= feedingPeriods[feedingPeriodPiv][1]) {
                        isFeeding = false;
                        feedingRateMult = 0;
                        feedingPeriodPiv++; // Update so we start looking for next period
                    } else
                        feedingRateMult = nominalFeedingRate;
                }

                if (includeHypoxiaAvoidance) {
                    checkAvoidanceCount++;
                    if (checkAvoidanceCount == checkAvoidanceInterval) {
                        checkAvoidanceCount = 0;
                        o2AffSum = SimUtils.setO2AffinityWithAvoidance(cageDims, dz, fishMaxDepth, mask, o2, affinity, o2Affinity);
                        //System.out.println("Updating O2 affinity. sum="+o2AffSum);

                    }
                }


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
                } else if (as.useTwin && (rank == N-1)) {
                    // This is the twin, introduce possible model error here.
                    for (int j = 0; j < ambientValueO2.length; j++) {
                        ambientValueO2_r[j] = ambientValueO2[j] + 0.25;

                    }
                    System.arraycopy(currentOffset, 0, currentOffset_r, 0, currentOffset.length);
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
                    SimUtils.interpolateVertical(interpAmbO2Par, as.parDepths, ambO2Par, cageDims[2], dz);

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


//                if (i==50) {
//                    double[] o2OutFlow2 = apOx.step(dt, o2, dxy, dz, useWalls, mask, 0, diffKappaO2, diffKappaO2Z,
//                            hydro, currentOffset_r, feedingRate, 0, ambientValueO2_r);
//
//                    System.out.println("i="+i+"\no2OutFlow vc: "+o2OutFlow[0]+"\no2OutFlow: "+o2OutFlow2[0]+"\n  Diff: "+(o2OutFlow[0]-o2OutFlow2[0]));
//                    System.out.println("1");
//                }

                double totalIntake = 0, o2ConsumptionRate = 0;
                double rho = 0;
                for (int ii=0; ii<cagePositions.size(); ii++) {
                    int[][] ranges = MultiCageUtils.getRanges(cagePositions.get(ii), dxy, rad);

                    double[] res = IngestionAndO2Subgrid.calculateIngestion(dt, fc, o2, affinity, o2Affinity,
                            o2AffSums[ii], ranges, ingDist, o2consDist, dxy, dz, mask, pelletWeight, ambientTemp, fish[ii], 0, useNewConsumptionModel);
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
                        ncfile.addGlobalAttribute("cageRad", spaghettiRad[0]);
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
                        SaveNetCDF.createProfileVariable(fishfile, "V", 0);
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

                    double[] groupArray = new double[fish[5].getNGroups()];

                    for (int j = 0; j < fishTmp.length; j++)
                        groupArray[j] = fish[0].getV(j);
                    SaveNetCDF.saveProfileVariable(fishfile, t, "V", 0, groupArray, true);
                    SaveNetCDF.saveScalarVariable(fishfile, t, "waste", outFlow, false);
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
                        groupArray[j] = fish[0].getAppetite(j);
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
                        currentComp[j] = currentProfile[j][0];
                    }
                    SaveNetCDF.saveProfileVariable(fishfile, t, "ext_currentU", 2, currentComp, false);
                    for (int j = 0; j < currentComp.length; j++) {
                        currentComp[j] = currentProfile[j][1];
                    }
                    SaveNetCDF.saveProfileVariable(fishfile, t, "ext_currentV", 2, currentComp, false);

                    // Save feeding rage (input value):
                    SaveNetCDF.saveScalarVariable(fishfile, t, "feedingRate", feedingRateMult, false);

                    // Save minimum O2 value:
                    //double[] values = SimUtils.minValueMeanAndFracHypoxia(o2, mask);
                    //SaveNetCDF.saveScalarVariable(fishfile, t, "min_O2", values[0], false);
                    //SaveNetCDF.saveScalarVariable(fishfile, t, "mean_O2", values[1], false);

                    //SaveNetCDF.saveScalarVariable(fishfile, t, "frac_hypoxia", values[2], false);

                    ArrayList<CageStats> cageStats = MultiCageUtils.getCageStats(o2, mask, cagePositions, rad, dxy,
                        SimUtils.HYPOXIA_THRESHOLD);
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

                    }


                    // Save o2 level at sensor positions:
                    for (int j=0; j<ms.names.length; j++) {
                        // Store value at sensor position:
                        SaveNetCDF.saveScalarVariable(fishfile, t, ms.names[j], o2[ms.pos[j][0]][ms.pos[j][1]][ms.pos[j][2]], false);
                        /*// Get values in a neighbourhood of the sensor to calculate spatial variability:
                        double nearStd = SimUtils.getStdAround(o2, o2Pos[j]);
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

}
