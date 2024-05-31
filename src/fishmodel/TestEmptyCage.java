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

import fishmodel.hydraulics.Balanced3DHydraulics;
import fishmodel.hydraulics.SimpleTankHydraulics;
import fishmodel.pellets.AdvectPellets;
import fishmodel.pellets.SimpleFish;
import fishmodel.hydraulics.CurrentMagicFields;
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
public class TestEmptyCage {

    public static final double HYPOXIA_THRESHOLD = 6;

    public static LinearInterpolator interpol = new LinearInterpolator();

    /**
     * Setup for Bjørøya
     */
    public static void main(String[] args) {

        // Save files:
        String saveDir = "./";
        String simNamePrefix = "empty_cm"; //"assim6_o2pert_lbeta_nopar_dropout";
        String simNamePostfix = "";


        // Simulation settings:
        boolean maskO2WhenSaving = false;
        boolean varyAmbient = false; // Reduction in ambient values towards the rest of the farm
        boolean decreasingCurrentFactor = true; // Model gradual decrease in current factor due to
                                                // increasing cage net biofouling
        boolean useCurrentMagic = true; // Use spatially variable current flow field
        CurrentMagicFields cmf = null;
        if (useCurrentMagic) {
            cmf = new CurrentMagicFields("C:/Users/alver/OneDrive - NTNU/prosjekt/O2_Bjørøya/currentmagic/currents_bjoroya_2m.nc");
        }

        boolean use3dBalancedCurrent = false; // Use Balanced3DHydraulics


        boolean useInstantaneousAmbientVals = true; // true to use ambient value time series, false to use daily averages


        // Simulation start time:
        int initYear = 2022, initMonth = Calendar.JUNE, initDate = 22, initHour = 0, initMin = 0, initSec = 0;
        double t_end = 2*3600;//1*24*3600; // Duration of simulation
        int nSim = 1; // Number of days to simulate (separate sims)
        int startAt = 0; // Set to >0 to skip one of more simulations, but count them in the sim numbering

        // Main current vector:
        double currU = 0.20, currV = 0.0;

        // Common settings:
        double rad = 25;
        double depth = 25, totDepth = 25; // Cage size (m)
        double dxy = 2, dz = dxy; // Model resolution (m)
        double dt = .5 * dxy; // Time step (s)
        int storeIntervalFeed = 600, storeIntervalInfo = 60;
        double fishMaxDepth = 20; // The maximum depth of the fish under non-feeding conditions

        double currentReductionFactor = ((useCurrentMagic || use3dBalancedCurrent) ? 1.0 : 0.8); // Multiplier for inside current
                                                                    // as function of outside

        // Environmental conditions:
        double currentSpeedInit = 2*0.04; // External current speed (m/s)
        double T_w = 14; double avO2 = 9; // mg / l
        //double T_w = 16; double avO2 = 8*0.9612; // mg / l
        //double T_w = 12; double avO2 = 8*1.0418; // mg / l

        // Oxygen diffusion constant (values updated further down)
        double diffKappaO2 = 0.1, diffKappaO2Z = 0.1;

        double[] currentOffset = new double[] {0,0,0}; // Global current vector

        // Fish setup (N, mean weight and std.dev weight):
        double nFishBjoroya = 0; // Estimated number of individuals in experimental period (source: FishTalk data)
        double meanWeight = 2869.5; // Estimated mean weight in experimental period (source: FishTalk data)
        double[] wFish = new double[] {meanWeight, 0.2*meanWeight};
        
        // Wind speed (x, y components in m/s) affecting feed spreader:
        double[] windSpeed = new double[] {0, 0};


        // Set up cage dimensions and cage grid:
        double modelDim = 2*(rad+4*dxy);
        //double modelDim = 2*(rad+2.5*rad); // TEST TEST TEST extra padding
        int[] cageDims = new int[3];
        cageDims[0] = (int)Math.ceil(modelDim/dxy);
        cageDims[1] = cageDims[0];
        cageDims[2] = (int)Math.ceil(depth/dz)+1;
        boolean[][][] mask = null;
        mask = CageMasking.circularMasking(cageDims, dxy, rad, false); // null

        boolean useWalls = false;
        System.out.println("Domain dimensions: ("+cageDims[0]+", "+cageDims[1]+", "+cageDims[2]+")");

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
        double nFish = nFishBjoroya;

        // Oxygen sensor positions:
        Measurements.MeasurementSet ms = Measurements.setupSensorPositionsBjoroya(cageDims, dxy, dz, rad);
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

        double[] affProfile = new double[] {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
        // Choose the flat or variable profile according to the "useVerticalDist" setting:
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
            ambientValueO2[i] = avO2;
        }

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
            // Here you set up the current profile (3D current vector per depth layer):
            double[][] currentProfile = new double[cageDims[2] + 1][3];
            for (int i=0; i<cageDims[2]+1; i++) {
                currentProfile[i][0] = 0;
            }
            hydro = SimpleTankHydraulics.getProfileHydraulicField(cageDims, currentProfile);

            if (decreasingCurrentFactor)
                currentReductionFactor = 0.8 + 0.05 - ((double)sim)*(0.2/*0.25*//8.0);


            AdvectPellets apOx = new AdvectPellets();
            if (varyAmbient) {
                apOx.setVaryAmbient(true, affinityProfile);

            }


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

            // Initialize O2 field:
            double extO2Val = 9;
            for (int i=0; i<cageDims[0]; i++)
                for (int j=0; j<cageDims[1]; j++)
                    for (int k=0; k<cageDims[2]; k++) {
                        o2[i][j][k] = extO2Val;
                    }


            // Initialize simple (grouped) fish model:
            SimpleFish fish = new SimpleFish(nFish, wFish[0], wFish[1]);
            double[][][] fishTmp = new double[fish.getNGroups()][1][1];

            SimpleDateFormat filenameForm = new SimpleDateFormat("dd_MM");
            String datePart = filenameForm.format(startTime);
            String filePrefix = simNamePrefix+datePart;

            // Establish file names to write data to:
            NetcdfFileWriteable ncfile = null;
            NetcdfFileWriteable fishfile = null;
            String ncfilePath = saveDir + filePrefix + simNamePostfix + ".nc";
            String fishfilePath = saveDir + filePrefix + simNamePostfix + "_fish.nc";
            boolean firstStore3d = true, firstStoreScalars = true;

            double totFeedAdded = 0;

            // Update variable current speed offset:
            currentOffset[0] = currU;//currentReductionFactor*currentSpeed*Math.cos(currentDirection*Math.PI/180.);
            currentOffset[1] = currV;//currentReductionFactor*currentSpeed*Math.sin(currentDirection*Math.PI/180.);
            double speed = Math.sqrt(currentOffset[0]*currentOffset[0] + currentOffset[1]*currentOffset[1]);
            double direc = Math.atan2(currentOffset[0], currentOffset[1])*180./Math.PI;
            double[] lDirections = new double[cageDims[2]],
                    lSpeeds = new double[cageDims[2]];
            for (int j = 0; j < cageDims[2]; j++) {
                lSpeeds[j] = speed;
                lDirections[j] = direc;

            }

            if (useCurrentMagic) {
                cmf.setCurrentField(hydro, lSpeeds, lDirections);
                currentOffset[0] = 0;
                currentOffset[1] = 0;
            }
            else if (use3dBalancedCurrent) {

                Balanced3DHydraulics.getTurbulentHydraulicField(cageDims, dxy, rad,
                        currentReductionFactor, currentProfile, hydro);
            }
            double t = 0;
            int n_steps = (int) (t_end / dt);
            long stime = System.currentTimeMillis();
            for (int i = 0; i < n_steps; i++) {

                double tMin = t / 60;




                diffKappaO2 = Math.min(0.5, 10*Math.pow(currentReductionFactor*0.06,2)); // Math.min(0.5, 10*Math.pow(currentReductionFactor*0.04,2));
                diffKappaO2Z = 5.0*0.1*Math.min(0.5, 10*Math.pow(currentReductionFactor*0.06,2)); // Math.min(0.5, 10*Math.pow(currentReductionFactor*0.04,2));
                //System.out.println("DiffKappa O2: "+diffKappaO2);

                setO2Boundaries(o2, extO2Val, extO2Val + 1.0);



                double[] o2OutFlow = apOx.step(dt, o2, dxy, dz, useWalls, mask, 0, diffKappaO2, diffKappaO2Z,
                        hydro, currentOffset, null, 0, ambientValueO2);




                double totalIntake = 0, rho = 0, o2ConsumptionRate = 0;

                t = t + dt;


                // Check if it is time to store 3D fields of feed and O2:
                if (i>0 && ((t/((double)storeIntervalFeed) - Math.floor(t/(double)storeIntervalFeed)) < 1e-5)) {
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
                if (i>0 && ((t/((double)storeIntervalInfo) - Math.floor(t/(double)storeIntervalInfo)) < 1e-5)) {

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
                    SaveNetCDF.saveScalarVariable(fishfile, t, "feedingRate", 0, false);

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

    public static void setO2Boundaries(double[][][] o2, double lowVal, double highVal) {
        int blocks = 6;

        // Left edge:
        for (int j=0; j<o2[0].length; j++) {
            for (int k = 0; k < o2[0][j].length; k++) {
                int jblock = j % blocks;
                boolean up = jblock == 0;
                o2[0][j][k] = up ? highVal : lowVal;
            }
        }

        // Lower edge:
        for (int i=0; i<o2.length; i++) {
            for (int k = 0; k < o2[i][0].length; k++) {
                int jblock = i % blocks;
                boolean up = jblock == 0;
                o2[i][0][k] = up ? highVal : lowVal;
            }
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
