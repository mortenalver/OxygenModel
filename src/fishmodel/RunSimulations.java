/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package fishmodel;

import fishmodel.hydraulics.SimpleTankHydraulics;
import fishmodel.pellets.*;

import java.io.File;
import java.io.IOException;
import java.text.NumberFormat;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.Date;
import java.util.Locale;
import java.util.Random;

import fishmodel.sim.CurrentMagicFields;
import fishmodel.sim.InputDataNetcdf;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import save.SaveNetCDF;
import ucar.nc2.NetcdfFileWriteable;

/**
 *
 * @author malv
 */
public class RunSimulations {

    public static final double HYPOXIA_THRESHOLD = 6;

    public static LinearInterpolator interpol = new LinearInterpolator();

    /**
     * Setup for Bjørøya
     */
    public static void main(String[] args) {

        NumberFormat nf1 = NumberFormat.getNumberInstance(Locale.ENGLISH);
        nf1.setMaximumFractionDigits(1);
        nf1.setMinimumFractionDigits(0);

        // Save files:
        String saveDir = "./";
        String simNamePrefix = "highdiff2_highero2_mag_2m_";
        String simNamePostfix = "";

        boolean[][][] mask = null;
        double lastMeanFeedDepth = -1;
        boolean maskO2WhenSaving = false;

        boolean varyAmbient = false; // Reduction in ambient values towards the rest of the farm
        double addRedMult = 0.015; // Scale factor for reduction in ambient values

        boolean useCurrentMagic = true;
        CurrentMagicFields cmf = null;
        if (useCurrentMagic) {
            cmf = new CurrentMagicFields("C:/Users/alver/OneDrive - NTNU/prosjekt/O2_Bjørøya/currentmagic/currents_bjoroya_2m.nc");
        }

        boolean useVerticalDist = false;

        boolean includeHypoxiaAvoidance = true;
        int checkAvoidanceInterval = 30, checkAvoidanceCount = 0;


        // Simulation start time:
        int initYear = 2022, initMonth = Calendar.JUNE, initDate = 22, initHour = 0, initMin = 0, initSec = 0;
        double t_end = 1*24*3600;//1*24*3600; // Duration of simulation
        int nSim = 9; // Number of days to simulate (separate sims)

        // Common settings:
        double rad = 25;
        int startAt = 0; // Set to >0 to skip one of more simulations, but count them in the sim numbering
        double depth = 25, totDepth = 25; // Cage size (m)
        double dxy = 2, dz = 2; // Model resolution (m)
        double dt = .5 * dxy; // Time step (s)
        int storeIntervalFeed = 600, storeIntervalInfo = 60;
        double fishMaxDepth = 20; // The maximum depth of the fish under non-feeding conditions

        double currentReductionFactor = (useCurrentMagic ? 1.0 : 0.8); // Multiplier for inside current as function of outside

        // Environmental conditions:
        double currentSpeedInit = 2*0.04; // External current speed (m/s)
        double T_w = 14; double avO2 = 9; // mg / l
        //double T_w = 16; double avO2 = 8*0.9612; // mg / l
        //double T_w = 12; double avO2 = 8*1.0418; // mg / l

        // Oxygen diffusion constant. To be set dependent on current speed.
        double diffKappaO2 = 0.1, diffKappaO2Z = 0.1;


        int[] cageDims = new int[3];
        double[] currentOffset = new double[] {0,0,0};

        // Oxygen sensor positions:
        String[] o2Names = new String[] {"C_5", "C_10", "C_15", "M1_5", "M1_10", "M1_15", "M2_5", "M2_10", "M2_15",
                "M3_5", "M3_10", "M3_15"};
        // angles: M1 128.2948, M2 2.8445, M3 246.8427
        double angle1 = 128.2948, angle2 = 2.8445, angle3 = 246.8427;
        double[] sensorAngles = new double[] {0, 0, 0, angle1, angle1, angle1, angle2, angle2, angle2,
                angle3, angle3, angle3};
        //double[] o2Rad = new double[] {0.6, 0.6, 0.6, 0.6}; // Distance from centre as fraction of cage radius
        double[] o2Rad = new double[] {0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1}; // Distance from centre as fraction of cage radius
        double[] o2Depth = new double[] {5, 10, 15, 5, 10, 15, 5, 10, 15, 5, 10, 15}; // Sensor depth (m)

        // Fish setup (N, mean weight and std.dev weight):
        double nFishBjoroya = 169821; // Estimated number of individuals in experimental period (source: FishTalk data)
        double meanWeight = 2869.5; // Estimated mean weight in experimental period (source: FishTalk data)
        double[] wFish = new double[] {meanWeight, 0.2*meanWeight};
        
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
        double modelDim = 2*(rad+4*dxy);
        cageDims[0] = (int)Math.ceil(modelDim/dxy);
        cageDims[1] = cageDims[0];
        cageDims[2] = (int)Math.ceil(depth/dz)+1;
        mask = CageMasking.circularMasking(cageDims, dxy, rad, false); // null
        boolean useWalls = false;

        System.out.println("Domain dimensions: ("+cageDims[0]+", "+cageDims[1]+", "+cageDims[2]+")");

        // Feeding setup:
        int[][] feedingPos = new int[][] {{cageDims[0]/2, cageDims[1]/2}};
        // Feeding periods (start/end in s):
        // Fra Eskil (Bjørøya): måltidene varte fra ca. kl. 07:30-17:30, i gjennomsnitt.
        int[][] feedingPeriods = new int[][] {{27000, 63000}, {86400+27000, 86400+63000}, {2*86400+27000, 2*86400+63000},
                {3*86400+27000, 3*86400+63000}, {4*86400+27000, 4*86400+63000}, {5*86400+27000, 5*86400+63000},
                {6*86400+27000, 6*86400+63000}, {7*86400+27000, 7*86400+63000}};
        int nPeriods = feedingPeriods.length;
        /*int[][] feedingPeriods = new int[][] {{1*3600, 2*3600}, {3*3600, 4*3600}, {5*3600, 6*3600},
                {7*3600, 8*3600}, {9*3600, 10*3600}, {11*3600, 12*3600}, {13*3600, 14*3600}, {15*3600, 16*3600},
                {17*3600, 18*3600}, {19*3600, 20*3600}, {21*3600, 22*3600}, {23*3600, 24*3600}, {25*3600, 26*3600},
                {27*3600, 28*3600}, {29*3600, 30*3600}, {31*3600, 32*3600}};*/
        for (int i=0; i<nPeriods; i++) {
            System.out.println("Feeding period "+(i+1)+": "+feedingPeriods[i][0]+" to "+feedingPeriods[i][1]);
        }
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

        // Determine number of fish, and feeding rate:
        double nFish = nFishBjoroya;
        System.out.println("N fish = "+nFish);
        double nominalFeedingRate = 2900.*1000/(10*3600); // Approximate feeding over 10 hours based on FishTalk data
        double feedingRateMult = 0; // Set each timestep
        System.out.println("Feeding rate = "+nominalFeedingRate);

        // Set up o2 sensor positions in grid:
        int[][] o2Pos = new int[o2Names.length][3];
        for (int i=0; i<o2Names.length; i++) {
            int xDist = (int)(o2Rad[i]*rad*Math.sin(Math.PI*sensorAngles[i]/180.)/dxy + cageDims[0]/2);
            int yDist = (int)(o2Rad[i]*rad*Math.cos(Math.PI*sensorAngles[i]/180.)/dxy + cageDims[1]/2);
            int zDist = (int)(o2Depth[i]/dz);
            System.out.println(o2Names[i]+": "+xDist+" , "+yDist+" , "+zDist);
            o2Pos[i][0] = xDist;
            o2Pos[i][1] = yDist;
            o2Pos[i][2] = zDist;
        }

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
        double[] affProfile_half = new double[] {0.5055, 0.5457, 0.9301, 1.5703, 1.8887, 1.8452, 1.7597, 1.6494, 1.5069,
                1.3724, 1.2942, 1.1834, 1.1174, 1.0362, 0.9690, 0.8882, 0.8552, 0.7947, 0.7804, 0.7334, 0.6966, 0.7004,
                0.6467, 0.5901, 0.5630, 0.5393, 0.5228, 0.5152};
        double[] affProfile_flat = new double[] {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
        double[] affProfile = useVerticalDist ? affProfile_half : affProfile_flat;


        double[] affDepths = new double[] {0.5000, 1.5000, 2.5000, 3.5000, 4.5000, 5.5000, 6.5000, 7.5000, 8.5000,
                9.5000, 10.5000, 11.5000, 12.5000, 13.5000, 14.5000, 15.5000, 16.5000, 17.5000, 18.5000, 19.5000,
                20.5000, 21.5000, 22.5000, 23.5000, 24.5000, 25.5000, 26.5000, 27.5000};
        /*double[] affProfile = new double[] {1, 1};
        double[] affDepths = new double[] {0, 30};*/
        double[] affinityProfile = new double[cageDims[2]];
        interpolateVertical(affinityProfile, affDepths, affProfile, cageDims[2], dz);

        for (int i = 0; i < affinityProfile.length; i++) {
            double v = affinityProfile[i];
            System.out.println("Affinity: "+v);
        }
        double[][][] o2Affinity = new double[cageDims[0]][cageDims[1]][cageDims[2]];
        double o2AffSum = setO2AffinityWithVerticalProfile(cageDims, dz, fishMaxDepth, mask, affinityProfile, o2Affinity, affinity);
        //double o2AffSum = setO2AffinityWithVerticalProfileAndEdgeDecrease(cageDims, rad, dxy, dz, fishMaxDepth, mask, affinityProfile, o2Affinity, affinity);
        int availableCellsForO2Uptake = countAvailableCellsForOxygenUptake(cageDims, dz, fishMaxDepth, mask);

        // Oxygen
        double[] ambientValueO2 = new double[cageDims[2]];
        for (int i = 0; i < ambientValueO2.length; i++) {
            ambientValueO2[i] = avO2;
        }

   
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
            c.add(Calendar.DAY_OF_MONTH, sim);
            Date startTime = c.getTime();

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
            if (varyAmbient) {
                apOx.setVaryAmbient(true, addRedMult,
                        new double[] {cageDims[0], cageDims[1]/2},
                        new double[] {1, 0});

            }

            // Initialize environmental input data:
            String inDataFile = "C:/Users/alver/OneDrive - NTNU/prosjekt/O2_Bjørøya/bjoroya_data.nc";
            if (!(new File(inDataFile)).exists())
                inDataFile = "bjoroya_data.nc";
            InputDataNetcdf inData = new InputDataNetcdf(inDataFile);
            inData.setStartTime(startTime);


            // Format a unit string for the time variable to save to NetCDF giving the initial time:
            SimpleDateFormat formatter = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
            String unitString = "seconds since "+formatter.format(startTime);
            System.out.println("Unit string: "+unitString);

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
            double[] ambVal = new double[] {inData.getO2Ambient5(), inData.getO2Ambient10(), inData.getO2Ambient15()};
            interpolateVertical(ambientValueO2, new double[] {5, 10, 15}, ambVal, cageDims[2], dz);
            for (int i=0; i<cageDims[0]; i++)
                for (int j=0; j<cageDims[1]; j++)
                    for (int k=0; k<cageDims[2]; k++) {
                        o2[i][j][k] = ambientValueO2[k];
                    }


            // Setup of surface feeding:
            double[][][] feedingRate = new double[cageDims[0]][cageDims[1]][cageDims[2]];
            double[][] fTemp = new double[cageDims[0]][cageDims[1]];

            double[][] surfFeed = new double[cageDims[0]][cageDims[1]];
            for (int i=0; i<feedingPos.length; i++) {
                PelletSpreaderModel.setPelletDist(fTemp, feedingPos[i][0], feedingPos[i][1], dxy, 0, 35, 0, true, windSpeed, null);
                // Add to total distribution:
                for (int ii=0; ii<surfFeed.length; ii++)
                    for (int jj = 0; jj < surfFeed[ii].length; jj++) {
                        surfFeed[ii][jj] = surfFeed[ii][jj] + fTemp[ii][jj];
                    }
            }

            for (int i=0; i<fTemp.length; i++)
                for (int j=0; j<fTemp[i].length; j++) {
                    feedingRate[i][j][0] = fTemp[i][j];
                }
            sourceTerm = feedingRate;

            // Initialize simple (grouped) fish model:
            SimpleFish fish = new SimpleFish(nFish, wFish[0], wFish[1]);
            double[][][] fishTmp = new double[fish.getNGroups()][1][1];

            SimpleDateFormat filenameForm = new SimpleDateFormat("dd_MM");
            String filePrefix = simNamePrefix+filenameForm.format(startTime);

            NetcdfFileWriteable ncfile = SaveNetCDF.initializeFile(saveDir + filePrefix + simNamePostfix+".nc", cageDims, 1, 1, unitString);
            NetcdfFileWriteable fishfile = SaveNetCDF.initializeFile(saveDir + filePrefix + simNamePostfix+"_fish.nc", new int[]{fish.getNGroups(), 1, cageDims[2]}, 1, 1, unitString);
            SaveNetCDF.createCageVariables(ncfile, "feed", "ingDist", "o2", "o2consDist");
            SaveNetCDF.createProfileVariable(fishfile, "appetite", 0);
            SaveNetCDF.createProfileVariable(fishfile, "ingested", 0);
            SaveNetCDF.createProfileVariable(fishfile, "V", 0);
            SaveNetCDF.createScalarVariables(fishfile, "rho", "feedingRate", "o2ConsumptionRate",
                    "min_O2", "mean_O2", "frac_hypoxia",
                    "meanFeedDepth", "d_meanFeedDepth", "totIngRate", "totIngested", "totFeed",
                    "waste", "waste_net");
            SaveNetCDF.createProfileVariable(fishfile, "ext_O2", 2); // dim=2 means along z dim
            SaveNetCDF.createProfileVariable(fishfile, "temperature", 2);
            SaveNetCDF.createProfileVariable(fishfile, "ext_currentU", 2);
            SaveNetCDF.createProfileVariable(fishfile, "ext_currentV", 2);
            SaveNetCDF.createScalarVariables(fishfile, o2Names);

            double totFeedAdded = 0;

            double currentDirection = 90;
            double currentSpeed = currentSpeedInit;
            Random rnd = new Random();
            double f_currentDir = Math.exp(-0.01*dt); // For simulating Gauss-Markov process
            double sigma_currentDir = 10;


            double t = 0;
            int n_steps = (int) (t_end / dt);
            long stime = System.currentTimeMillis();
            for (int i = 0; i < n_steps; i++) {

                //System.out.println("t = "+t);
                double tMin = t / 60;

                if (inData.advance(t) || (i==0)) {
                    //System.out.println("Updating environment: t = "+t);

                    double[] tempVal = {inData.getTemperature5(), inData.getTemperature10(), inData.getTemperature15()};
                    //tempVal[0] = 18; tempVal[1] = 10; tempVal[2] = 2;
                    interpolateVertical(ambientTemp, new double[] {5, 10, 15}, tempVal, cageDims[2], dz);
                    ambVal = new double[] {inData.getO2Ambient5(), inData.getO2Ambient10(), inData.getO2Ambient15()};
                    //ambVal[0] = 6; ambVal[1] = 7; ambVal[2] = 8;
                    //double[] ambVal = {2., 10., 5.};
                    interpolateVertical(ambientValueO2, new double[] {5, 10, 15}, ambVal, cageDims[2], dz);
                    /*for (int j = 0; j < tempVal.length; j++) {
                        double v = tempVal[j];
                        System.out.println("Temp val: "+v);
                    }
                    for (int j = 0; j < ambientTemp.length; j++) {
                        double v = ambientTemp[j];
                        System.out.println((i+1)+": "+v);
                    }*/


                    //currentSpeed = inData.getExtCurrentSpeed();
                    //currentDirection = inData.getExtCurrentDir();
                    double[] obsCurrentProfile = inData.getExtCurrentSpeedProfile();
                    double[] obsCurrentDirProfile = inData.getExtCurrentDirProfile();
                    double[] obsCurrentComp1 = new  double[obsCurrentProfile.length],
                            obsCurrentComp2 = new  double[obsCurrentProfile.length];
                    // Current directions are given as the direction the current flows towards, with
                    // 0 degrees being north and 90 degrees being east. Verified by comparing histograms
                    // with the textual descriptions in the report by Aqua Kompetanse.
                    // x component: speed*sin(direction)
                    // y component: speed*cos(direction)
                    for (int j = 0; j < obsCurrentComp1.length; j++) {
                        obsCurrentComp1[j] = currentReductionFactor*
                                obsCurrentProfile[j]*Math.sin(obsCurrentDirProfile[j]*Math.PI/180.);
                        obsCurrentComp2[j] = currentReductionFactor*
                                obsCurrentProfile[j]*Math.cos(obsCurrentDirProfile[j]*Math.PI/180.);
                    }

                    double[] obsCurrentDepths = inData.getCurrentDepths();
                    double[] interpProfile1 = new double[cageDims[2]],
                            interpProfile2 = new double[cageDims[2]];
                    interpolateVertical(interpProfile1, obsCurrentDepths, obsCurrentComp1, cageDims[2], dz);
                    interpolateVertical(interpProfile2, obsCurrentDepths, obsCurrentComp2, cageDims[2], dz);

                    if (!useCurrentMagic) {
                        for (int j = 0; j < interpProfile1.length; j++) {
                            //System.out.println("Interpolated current speed "+(j)+": "+interpProfile1[j]+" / "+interpProfile2[j]);
                            currentProfile[j][0] = interpProfile1[j];
                            currentProfile[j][1] = interpProfile2[j];
                            currentProfile[j][2] = 0.;
                        }
                        // Update current field using the new profile:
                        SimpleTankHydraulics.getProfileHydraulicField(hydro, cageDims, currentProfile);
                    }
                    else {
                        double[] lDirections = new double[cageDims[2]],
                                lSpeeds = new double[cageDims[2]];
                        for (int j = 0; j < cageDims[2]; j++) {
                            lSpeeds[j] = Math.sqrt(interpProfile1[j]*interpProfile1[j] + interpProfile2[j]*interpProfile2[j]);
                            lDirections[j] = Math.atan2(interpProfile1[j], interpProfile2[j])*180./Math.PI;
                            //System.out.println(nf.format(interpProfile1[j])+" / "+nf.format(interpProfile2[j])+" --> speed="+
                            //        nf.format(lSpeeds[j])+", dir="+nf.format(lDirections[j]));
                        }
                        cmf.setCurrentField(hydro, lSpeeds, lDirections);

                    }
                }

                // Update variable current speed offset:
                currentOffset[0] = 0;//currentReductionFactor*currentSpeed*Math.cos(currentDirection*Math.PI/180.);
                currentOffset[1] = 0;//currentReductionFactor*currentSpeed*Math.sin(currentDirection*Math.PI/180.);

                diffKappaO2 = Math.min(0.5, 10*Math.pow(currentReductionFactor*0.06,2)); // Math.min(0.5, 10*Math.pow(currentReductionFactor*0.04,2));
                diffKappaO2Z = 5.0*0.1*diffKappaO2;
                //System.out.println("DiffKappa O2: "+diffKappaO2);

                // Update feeding rate depending on preset feeding periods:
                if (!isFeeding) { // Not already feeding. Check if we should start:
                    if ((feedingPeriodPiv < feedingPeriods.length) && (t >= feedingPeriods[feedingPeriodPiv][0])) {
                        isFeeding = true;
                        feedingRateMult = nominalFeedingRate;
                        // Reset gut content at start of feeding period:
                        fish.resetAllV();
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

                /*if (includeHypoxiaAvoidance) {
                    checkAvoidanceCount++;
                    if (checkAvoidanceCount == checkAvoidanceInterval) {
                        checkAvoidanceCount = 0;
                        o2AffSum = setO2AffinityWithAvoidance(cageDims, dz, fishMaxDepth, mask, o2, affinity, o2Affinity);
                        //System.out.println("Updating O2 affinity. sum="+o2AffSum);

                    }
                }*/


                totFeedAdded += dt * feedingRateMult;

                double[] r = ap.step(dt, fc, dxy, dz, useWalls, mask, sinkingSpeed, diffKappa, diffKappaZ, 
                        hydro, currentOffset, sourceTerm, feedingRateMult, ambientValueFeed);
                outFlow = r[0]; // Feed lost from grid (not used)
                outFlow_net = r[1]; // Feed lost from the unmasked part of the grid (feed lost through side)

                double[] o2OutFlow = apOx.step(dt, o2, dxy, dz, useWalls, mask, 0, diffKappaO2, diffKappaO2Z,
                        hydro, currentOffset, feedingRate, 0, ambientValueO2);

                double[] res = IngestionAndO2Tempprofile.calculateIngestion(dt, fc, o2, affinity, o2Affinity, o2AffSum,
                        availableCellsForO2Uptake, ingDist, o2consDist, dxy, dz, mask, pelletWeight, ambientTemp, fish);
                double totalIntake = res[0], rho = res[1], o2ConsumptionRate = res[2];

                t = t + dt;

                if (i>0 && ((t/((double)storeIntervalFeed) - Math.floor(t/(double)storeIntervalFeed)) < 1e-5)) {
                    double elapsed = (double) ((System.currentTimeMillis() - stime)) / 60000.;
                    double fractionCompleted = ((double) i) / ((double) n_steps);
                    double remaining = (elapsed / fractionCompleted) - elapsed;
                    System.out.println("t = " + nf1.format(t) + " - Estimated time to complete: " + nf1.format(remaining) + " minutes");

                    SaveNetCDF.saveCageVariable(ncfile, t, "feed", fc, mask, true);
                    SaveNetCDF.saveCageVariable(ncfile, t, "ingDist", ingDist, mask, false);
                    SaveNetCDF.saveCageVariable(ncfile, t, "o2", o2, (maskO2WhenSaving ? mask : null), false);
                    SaveNetCDF.saveCageVariable(ncfile, t, "o2consDist", o2consDist, mask, false);
                }

                if (i>0 && ((t/((double)storeIntervalInfo) - Math.floor(t/(double)storeIntervalInfo)) < 1e-5)) {
                //if ((t - Math.floor(t) < dt) && (Math.floor(t) % storeIntervalInfo == 0)) {

                    double[] groupArray = new double[fish.getNGroups()];

                    for (int j = 0; j < fishTmp.length; j++)
                        groupArray[j] = fish.getV(j);
                    SaveNetCDF.saveProfileVariable(fishfile, t, "V", 0, groupArray, true);
                    SaveNetCDF.saveScalarVariable(fishfile, t, "waste", outFlow, false);
                    SaveNetCDF.saveScalarVariable(fishfile, t, "waste_net", outFlow_net, false);
                    
                    double totI = 0, totIngRate = 0;
                    for (int j = 0; j < fishTmp.length; j++) {
                        groupArray[j] = fish.getIngested(j);
                        totI += fish.getN(j) * fishTmp[j][0][0];
                        totIngRate += fish.getN(j) * fish.getIngRate(j);
                    }
                    SaveNetCDF.saveProfileVariable(fishfile, t, "ingested", 0, groupArray, false);


                    for (int j = 0; j < groupArray.length; j++) {
                        groupArray[j] = fish.getAppetite(j);
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
                    double[] values = minValueMeanAndFracHypoxia(o2, mask);
                    SaveNetCDF.saveScalarVariable(fishfile, t, "min_O2", values[0], false);
                    SaveNetCDF.saveScalarVariable(fishfile, t, "mean_O2", values[1], false);
                    SaveNetCDF.saveScalarVariable(fishfile, t, "frac_hypoxia", values[2], false);

                    // Save o2 level at sensor positions:
                    for (int j=0; j<o2Names.length; j++) {
                        // Store value at sensor position:
                        SaveNetCDF.saveScalarVariable(fishfile, t, o2Names[j], o2[o2Pos[j][0]][o2Pos[j][1]][o2Pos[j][2]], false);
                        /*// Get values in a neighbourhood of the sensor to calculate spatial variability:
                        double nearStd = getStdAround(o2, o2Pos[j]);
                        SaveNetCDF.saveScalarVariable(fishfile, t, o2Names[j]+"_std", nearStd, false);*/
                    }

                    lastMeanFeedDepth = meanFeedDepth;

                }


            }


            double totI = 0;
            for (int j = 0; j < fishTmp.length; j++) {
                totI += fish.getN(j) * fish.getIngested(j);
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
