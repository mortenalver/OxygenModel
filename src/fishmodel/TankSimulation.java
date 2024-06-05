/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package fishmodel;


import fishmodel.hydraulics.StationaryTankFlow;
import fishmodel.pellets.AdvectPellets;
import fishmodel.pellets.IngestionAndO2Tempprofile;
import fishmodel.pellets.SimpleFish;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import save.SaveNetCDF;
import ucar.nc2.NetcdfFileWriteable;

import java.io.IOException;
import java.text.NumberFormat;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.Date;
import java.util.Locale;
import java.util.Random;


/**
 *
 * @author malv
 */
public class TankSimulation {

    public static final double HYPOXIA_THRESHOLD = 6;

    public static LinearInterpolator interpol = new LinearInterpolator();

    /**
     * Setup for Bjørøya
     */
    public static void main(String[] args) {

        // Save files:
        String saveDir = "./";
        String simNamePrefix = "tanktest4";
        String simNamePostfix = "";


        // Simulation settings:
        boolean maskO2WhenSaving = false;

        boolean useVerticalDist = false; // Use non-uniform vertical distribution (defined further down)
                                        // for non-feeding fish


        //boolean includeHypoxiaAvoidance = true;     int checkAvoidanceInterval = 30, checkAvoidanceCount = 0;

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

        // Simulation start time:
        int initYear = 2022, initMonth = Calendar.JUNE, initDate = 27, initHour = 0, initMin = 0, initSec = 0;
        initDate += daysToAdd;
        double t_end = 60*60+1;//1*24*3600; // Duration of simulation
        int nSim = 1; // Number of days to simulate (separate sims)
        int startAt = 0; // Set to >0 to skip one of more simulations, but count them in the sim numbering

        // Common settings:
        double rad = 3;
        double depth = 1.1, totDepth = 1; // Cage size (m)
        double dxy = 0.1, dz = dxy; // Model resolution (m)
        double dt = .5 * dxy; // Time step (s)
        double storeIntervalSFeed = 60, storeIntervalSInfo = 60;
        int storeIntervalFeed = (int)Math.round(storeIntervalSFeed/dt);
        int storeIntervalInfo = (int)Math.round(storeIntervalSInfo/dt);
        System.out.println("Store intervals: "+storeIntervalFeed+" / "+storeIntervalInfo);
        double fishMaxDepth = 20; // The maximum depth of the fish under non-feeding conditions

        double currentReductionFactor = 1.0; // Multiplier for inside current as function of outside

        // Environmental conditions:
        double T_w = 14; double avO2 = 8; // mg / l

        // Oxygen diffusion constant (values updated further down)
        double diffKappaO2 = 0.1, diffKappaO2Z = 0.1;

        double[] currentOffset = new double[] {0,0,0}; // Global current vector
        double[] currentOffset_r = new double[] {0,0,0}; // Perturbed global current vector

        // Fish setup (N, mean weight and std.dev weight):
        double nFishTank = 1e-7*150; // Estimated number of individuals in experimental period (source: FishTalk data)
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
        double modelDim = 2*(rad+6*dxy);
        //double modelDim = 2*(rad+2.5*rad); // TEST TEST TEST extra padding
        int[] cageDims = new int[3];
        cageDims[0] = (int)Math.ceil(modelDim/dxy);
        cageDims[1] = cageDims[0];
        cageDims[2] = (int)Math.ceil(depth/dz)+1;
        boolean[][][] mask = null;
        mask = CageMasking.circularMasking(cageDims, dxy, rad, true); // null

        boolean useWalls = true;
        System.out.println("Domain dimensions: ("+cageDims[0]+", "+cageDims[1]+", "+cageDims[2]+")");


        // Feeding setup:
        // Feeding periods (start/end in s):
        // Fra Eskil (Bjørøya): måltidene varte fra ca. kl. 07:30-17:30, i gjennomsnitt.
        int[][] feedingPeriods = new int[][] {{0, 3600}};
        int nPeriods = feedingPeriods.length;
        for (int i=0; i<nPeriods; i++) {
            System.out.println("Feeding period "+(i+1)+": "+feedingPeriods[i][0]+" to "+feedingPeriods[i][1]);
        }
        Object sourceTerm = null;

        double[] ambientValueFeed = new double[cageDims[2]]; // Outside feed concentrations are set to 0
        for (int i = 0; i < ambientValueFeed.length; i++) {
            ambientValueFeed[i] = 0;
        }

        // Set up array for ambient temperature:
        double[] ambientTemp = new double[cageDims[2]];
        for (int i = 0; i < ambientTemp.length; i++) {
            ambientTemp[i] = T_w;
        }

        // Determine number of fish, and feeding rate:
        double nFish = nFishTank;
        System.out.println("N fish = "+nFish);
        double nominalFeedingRate = 2900.*1000/(10*3600); // Approximate feeding over 10 hours based on FishTalk data
        double feedingRateMult = 0; // Set each timestep
        System.out.println("Feeding rate = "+nominalFeedingRate);

        // Oxygen sensor positions:
        Measurements.MeasurementSet ms = Measurements.setupSensorPositionsTank(cageDims, dxy, dz, rad);
        String[] o2Names = ms.names;
        int[][] o2Pos = ms.pos;
        // List O2 sensor positions and which are masked / not masked:
        for (int i=0; i<o2Names.length; i++) {
            System.out.println(o2Names[i]+": "+o2Pos[i][0]+" , "+o2Pos[i][1]+" , "+o2Pos[i][2]+", mask="+mask[o2Pos[i][0]][o2Pos[i][1]][o2Pos[i][2]]);
        }

        // Feed affinity (determines how "interested" the fish is in feeding in each model grid cell):
        double[][][] affinity = new double[cageDims[0]][cageDims[1]][cageDims[2]];
        // O2 affinity (determines the typical vertical distribution of the fish):
        double[][][] o2Affinity = new double[cageDims[0]][cageDims[1]][cageDims[2]];

        // Define vertical distribution based on telemetry data (from 8 individuals):
        double[] affProfile_orig = new double[] {0.0110, 0.0913, 0.8601, 2.1406, 2.7774, 2.6903, 2.5195, 2.2987, 2.0137,
                1.7448, 1.5883, 1.3667, 1.2348, 1.0724, 0.9379, 0.7764, 0.7104, 0.5895, 0.5607, 0.4668, 0.3933,
                0.4009, 0.2935, 0.1801, 0.1260, 0.0787, 0.0457, 0.0304};
        double[] affProfile_flat = new double[] {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
        // Choose the flat or variable profile according to the "useVerticalDist" setting:
        double[] affProfile = useVerticalDist ? affProfile_orig : affProfile_flat;
        // Depths values for the values defined above:
        double[] affDepths = new double[] {0.5000, 1.5000, 2.5000, 3.5000, 4.5000, 5.5000, 6.5000, 7.5000, 8.5000,
                9.5000, 10.5000, 11.5000, 12.5000, 13.5000, 14.5000, 15.5000, 16.5000, 17.5000, 18.5000, 19.5000,
                20.5000, 21.5000, 22.5000, 23.5000, 24.5000, 25.5000, 26.5000, 27.5000};

        double[] affinityProfile = new double[cageDims[2]];
        interpolateVertical(affinityProfile, affDepths, affProfile, cageDims[2], dz);

        // Set feeding and O2 affinity based on the chosen vertical profile and cage mask:
        double o2AffSum = setO2AffinityWithVerticalProfile(cageDims, dz, fishMaxDepth, mask, affinityProfile, o2Affinity, affinity);
        //double o2AffSum = setO2AffinityWithVerticalProfileAndEdgeDecrease(cageDims, rad, dxy, dz, fishMaxDepth, mask, affinityProfile, o2Affinity, affinity);
        int availableCellsForO2Uptake = countAvailableCellsForOxygenUptake(cageDims, dz, fishMaxDepth, mask);

        // Oxygen
        double[] ambientValueO2 = new double[cageDims[2]]; // Ambient value of O2
        for (int i = 0; i < ambientValueO2.length; i++) {
            ambientValueO2[i] = 20;///avO2;
        }

        avO2 = 10;

        // Initialize number formatter:
        NumberFormat nf1 = NumberFormat.getNumberInstance(Locale.ENGLISH);
        nf1.setMaximumFractionDigits(1); nf1.setMinimumFractionDigits(0);
        // Initialize random number generator:
        Random rnd = new Random();
        // Other initialization:
        double lastMeanFeedDepth = -1;


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
            StationaryTankFlow tankFlow = new StationaryTankFlow("C:\\Users\\alver\\OneDrive - NTNU\\prosjekt\\AQUAEXCEL3\\particleTransport\\inrae_orig_max_3d_4.nc",
                    cageDims, dxy, dz);
            hydro = tankFlow.getFlowField();
            //int[] inletPos = tankFlow.getCellForCoordinates(new double[] {3, 0.15, 0.77}, dxy, dz);
            //int[] inletPos = tankFlow.getCellForCoordinates(new double[] {4.16, 1.56, 0.142}, dxy, dz);
            int[] inletPos = new int[] {7, 36, 2};
            System.out.println("Inlet cell: "+inletPos[0]+", "+inletPos[1]+", "+inletPos[2]+", mask="+mask[inletPos[0]][inletPos[1]][inletPos[2]]);
            System.out.println("w curr at inlet: "+hydro[inletPos[0]][inletPos[1]][inletPos[2]][2]+
                    " , "+hydro[inletPos[0]][inletPos[1]][inletPos[2]+1][2]);

            AdvectPellets ap = new AdvectPellets();
            AdvectPellets apOx = new AdvectPellets();


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
            for (int i=0; i<cageDims[0]; i++)
                for (int j=0; j<cageDims[1]; j++)
                    for (int k=0; k<cageDims[2]; k++) {
                        o2[i][j][k] = avO2;
                    }

            // Setup of surface feeding:
            double[][][] feedingRate = new double[cageDims[0]][cageDims[1]][cageDims[2]];
            double[][] fTemp = new double[cageDims[0]][cageDims[1]];

            double[][] surfFeed = new double[cageDims[0]][cageDims[1]];
            int nCells = 0;
            int centerX = (int)(cageDims[0]/3), centerY = cageDims[1]/2;
            for (int i=0; i<surfFeed.length; i++)
                for (int j=0; j<surfFeed[i].length; j++) {
                    if (mask[i][j][0]) {
                        double dist = Math.sqrt((double)(i-centerX)*(double)(i-centerX)
                                + (double)(j-centerY)*(double)(j-centerY));
                        if (dist < 5) {
                            surfFeed[i][j] = 1;
                            nCells++;
                        }
                    }
                }
            double divideBy = (double)nCells;
            for (int i=0; i<surfFeed.length; i++)
                for (int j=0; j<surfFeed[i].length; j++) {
                    if (mask[i][j][0]) {
                        feedingRate[i][j][0] = surfFeed[i][j]/divideBy;
                    }
                }

            sourceTerm = feedingRate;

            // Initialize simple (grouped) fish model:
            SimpleFish fish = new SimpleFish(nFish, wFish[0], wFish[1]);
            double[][][] fishTmp = new double[fish.getNGroups()][1][1];

            SimpleDateFormat filenameForm = new SimpleDateFormat("dd_MM");
            String datePart = filenameForm.format(startTime);
            String filePrefix = simNamePrefix+datePart;

            // Establish file names to write data to:
            NetcdfFileWriteable ncfile = null;
            NetcdfFileWriteable fishfile = null;
            String ncfilePath = saveDir + filePrefix + simNamePostfix +".nc";
            String fishfilePath = saveDir + filePrefix + simNamePostfix +"_fish.nc";
            boolean firstStore3d = true, firstStoreScalars = true;


            double totFeedAdded = 0;


            double t = 0;
            int n_steps = (int) (t_end / dt);
            long stime = System.currentTimeMillis();
            for (int i = 0; i < n_steps; i++) {

                //System.out.println("t = "+t);
                double tMin = t / 60;


                // Update variable current speed offset:
                currentOffset[0] = 0.;//currentReductionFactor*currentSpeed*Math.cos(currentDirection*Math.PI/180.);
                currentOffset[1] = 0.;//currentReductionFactor*currentSpeed*Math.sin(currentDirection*Math.PI/180.);

                diffKappaO2 = 0.001;
                diffKappaO2Z = diffKappaO2;
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


                double maxval = 0;


                double[] r = ap.step(dt, fc, dxy, dz, useWalls, mask, sinkingSpeed, diffKappa, diffKappaZ,
                        hydro, currentOffset_r, sourceTerm, feedingRateMult, ambientValueFeed);
                outFlow = r[0]; // Feed lost from grid (not used)
                outFlow_net = r[1]; // Feed lost from the unmasked part of the grid (feed lost through side)


                // Set inlet o2 value at inlet position each time step:
                o2[inletPos[0]][inletPos[1]][inletPos[2]] += dt*0.5;

                double[] o2OutFlow = apOx.step(dt, o2, dxy, dz, useWalls, mask, 0, diffKappaO2, diffKappaO2Z,
                        hydro, currentOffset_r, feedingRate, 0, ambientValueO2);


                double[] res = IngestionAndO2Tempprofile.calculateIngestion(dt, fc, o2, affinity, o2Affinity, o2AffSum,
                        availableCellsForO2Uptake, ingDist, o2consDist, dxy, dz, mask, pelletWeight, ambientTemp, fish,
                        0);
                double totalIntake = res[0], rho = res[1], o2ConsumptionRate = res[2];



                t = t + dt;




                // Check if it is time to store 3D fields of feed and O2:
                if (i>0 && ((i+1)%storeIntervalFeed == 0)) {
                //((t/((double)storeIntervalFeed) - Math.floor(t/(double)storeIntervalFeed)) < 1e-5)) {
                    double elapsed = (double) ((System.currentTimeMillis() - stime)) / 60000.;
                    double fractionCompleted = ((double) i) / ((double) n_steps);
                    double remaining = (elapsed / fractionCompleted) - elapsed;
                    System.out.println("t = " + nf1.format(t) + " - Estimated time to complete: " + nf1.format(remaining) + " minutes");

                    if (firstStore3d) {
                        firstStore3d = false;
                        ncfile = SaveNetCDF.initializeFile(ncfilePath, cageDims, 1, 1, unitString, ms);
                        SaveNetCDF.createCageVariables(ncfile, "feed", "ingDist", "o2", "o2consDist");
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

                // Check if it is time to store scalar output values:
                if (i>0 && ((i+1)%storeIntervalInfo == 0)) {
                    //(i>0 && ((t/((double)storeIntervalInfo) - Math.floor(t/(double)storeIntervalInfo)) < 1e-5)) {

                    if (firstStoreScalars) {
                        firstStoreScalars = false;
                        fishfile = SaveNetCDF.initializeFile(fishfilePath, new int[]{fish.getNGroups(), 1, cageDims[2]}, 1, 1, unitString, ms);
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
                    }
                    else {
                        try {
                            fishfile = NetcdfFileWriteable.openExisting(fishfilePath);
                        } catch (IOException e) {
                            e.printStackTrace();
                        }
                    }

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

                    try {
                        fishfile.close();
                    } catch (IOException e) {
                        e.printStackTrace();
                    }

                }


            }


            double totI = 0;
            for (int j = 0; j < fishTmp.length; j++) {
                totI += fish.getN(j) * fish.getIngested(j);
            }

            System.out.println("totFeedAdded = " + totFeedAdded);
            System.out.println("totI = " + totI);
            System.out.println("Feed wastage: " + nf.format(100 * (totFeedAdded - totI) / totFeedAdded) + " %");

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

    private static double setO2AffinityWithVerticalProfileAndDirection(int[] cageDims, double dz, double fishMaxDepth, double[] dirVector, boolean[][][] mask,
                                                           double[] affinityProfile, double[][][] affinity,
                                                           double[][][] o2Affinity) {
        double o2AffSum = 0;
        double centerX = ((double)cageDims[0])/2.0,
                centerY = ((double)cageDims[1])/2.0;
        double dirVectorLength = Math.sqrt(dirVector[0]*dirVector[0] + dirVector[1]*dirVector[1]);

        for (int i=0; i<cageDims[0]; i++)
            for (int j=0; j<cageDims[1]; j++) {
                // Set up unit vector in the direction from center to this horizontal position:
                double distX = ((double)i) - centerX,
                        distY = ((double)j) - centerY;
                double distTot = Math.sqrt(distX*distX+distY*distY);
                if (distTot > 0) {
                    distX = distX/distTot;
                    distY = distY/distTot;
                }
                // Calculate dot product between direction vector and the input dir vector. This gives a value
                // that is equal to the length of the input dir vector multiplied by the cosine of the angle between
                // the two vectors:
                double dotProduct = distX*dirVector[0] + distY*dirVector[1];
                //dotProduct = (dotProduct+1)*(dotProduct+1)-1;
                System.out.println("distX="+distX+", distY="+distY+", dotProd="+dotProduct);


                for (int k = 0; k < cageDims[2]; k++) {
                    double lDepth = ((double) k + 0.5) * dz;
                    if ((mask == null || mask[i][j][k]) && (lDepth < fishMaxDepth)) {
                        // Check oxygen avoidance criterion:
                        //affinity[i][j][k] = Math.max(0., affinityProfile[k] + dotProduct);
                        //o2Affinity[i][j][k] = Math.max(0., affinityProfile[k] + dotProduct);
                        if (dotProduct/dirVectorLength > 0*0.2588) {
                            affinity[i][j][k] = affinityProfile[k];
                            o2Affinity[i][j][k] = affinityProfile[k];
                        } else {
                            affinity[i][j][k] = affinityProfile[k];
                            o2Affinity[i][j][k] = 0;
                        }
                        //affinity[i][j][k] = affinityProfile[k];
                        //o2Affinity[i][j][k] = affinityProfile[k];

                    } else {
                        affinity[i][j][k] = 0;
                        o2Affinity[i][j][k] = 0;

                    }

                    o2AffSum += o2Affinity[i][j][k];
                }
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
