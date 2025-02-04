package fishmodel.pellets;

import fishmodel.CageMasking;
import fishmodel.hydraulics.SimpleTankHydraulics;
import save.SaveNetCDF;
import ucar.ma2.Array;
import ucar.nc2.NetcdfFileWriteable;

import java.io.IOException;
import java.text.NumberFormat;
import java.util.Arrays;
import java.util.Locale;


public class AdvectPellets {

    static int numProcessors = Math.min(40, Runtime.getRuntime().availableProcessors());

    static int FEEDING_KLEVEL = 0;

    int counter = 0;

    private boolean varyAmbient = false;
    //private double[] ambRedCenter = null;
    //private double[] dsVector = null;
    //private double ambRedMult = 0;
    private double[] affinityProfile;

    private double[][][] advect = null, diffus = null, newValues = null;

    public static void disableMultiprocessing() {
        numProcessors = 1;
    }

    /**
     * Steps pellet model.
     * @param dt The time step.
     * @param pConc The pellet distribution to be updated.
     * @param dxy The horizontal grid cell width.
     * @param dz The vertical grid cell depth.
     * @param useWalls If true, and a cage mask is given, simulate a setting with walls determined by the mask.
     *    If false, and a cage mask is given, simulate a cage setting where the mask has no effect on transport,
     *    but feed passing outside the cage is considered wasted in the calculation of feed loss.
     * @param mask Mask denoting which cells are inside a wall. A null value signals that no walls are present.
     * @param sinkingSpeed The sinking speed of pellets.
     * @param diffKappaXY The horizontal diffusion parameter kappa.
     * @param diffKappaZ The vertical diffusion parameter kappa. Note that for feed pellets it makes sense to use a
     *                   significantly larger value for vertical diffusion than for horizontal, in order to get
     *                   results harmonizing with pellet drop and horizontal spread experiments performed by
     *                   Kristoffer Skøien et al. - this is likely to be primarily caused by the spread in
     *                   sinking speed of pellets rather than diffusion/turbulence.
     *                   The initial version of the model suffered from numerical diffusion, and in the calibration
     *                   simulations no horizontal speed was used, thereby cancelling numerical diffusion in the
     *                   horizontal plane. This lead to the erroneous conclusion that a single kappa value gave
     *                   good match both horizontally and vertically.
     *                   After implementing an advection scheme suppressing numerical diffusion, the calibration
     *                   runs were repeated, and a factor of 25 greater vertical kappa was found to give suitable
     *                   match with the experimental data.
     * @param currentSpeed The (3d) current speed in all cells.
     * @param currentSpeedOffset Offset added to current speed
     * @param sourceTerm The feeding rate - either 2D array of feeding rate distributed over all surface cells,
     *                   or 3D array of feeding rate distributed over whole grid
     * @param feedingRateMult A multiplier for the feeding rate.
     * @param ambientValue The depth-dependent outside value to apply if edges are open.
     * @return The amount of feed lost through the outer grid edges.
     */
    public double[] step(final double dt, final double[][][] pConc, final double dxy, final double dz, final boolean useWalls,
                       final boolean[][][] mask, final double sinkingSpeed, final double diffKappaXY, final double diffKappaZ, 
                       final double[][][][] currentSpeed, double[] currentSpeedOffset, 
                       final Object sourceTerm, final double feedingRateMult, final double[] ambientValue) {

        counter++;
        // Find out what sort of feed addition rate we have:
        final double[][] feedingRate;
        final double[][][] distFeeding;
        if (sourceTerm instanceof double[][])
            feedingRate = (double[][])sourceTerm;
        else feedingRate = null;
        if (sourceTerm instanceof double[][][])
            distFeeding = (double[][][])sourceTerm;
        else distFeeding = null;

        // Calculate advection and diffusion, taking the mask into account only if we should simulate a tank:
        final int[] dim = new int[] {pConc.length, pConc[0].length, pConc[0][0].length};
        if (advect == null) {
            advect = new double[dim[0]][dim[1]][dim[2]];
            diffus = new double[dim[0]][dim[1]][dim[2]];
            newValues = new double[dim[0]][dim[1]][dim[2]];
        }

        //final double[][][] presu = new double[dim[0]][dim[1]][dim[2]];
        //final double[][][] posu = new double[dim[0]][dim[1]][dim[2]];

        class CalcPart implements Runnable {
            int kstart, kend;
            public double preSum = 0, postSum = 0, preNot = 0, postNot = 0;
            public CalcPart(int kstart, int kend) {
                this.kstart = kstart;
                this.kend = kend;
                //System.out.println("kstart = "+kstart+", kend = "+kend);
            }
            @Override
            public void run() {

                //calcAdvectAndDiff(dim, kstart, kend, pConc, dxy, dz, useWalls ? mask : null, sinkingSpeed, diffKappaXY, currentSpeed,
                //        currentSpeedOffset, ambientValue);
                calcAdvectAndDiff2(dim, kstart, kend, pConc, dxy, dz, useWalls ? mask : null, sinkingSpeed, diffKappaXY, diffKappaZ,
                        currentSpeed, currentSpeedOffset, ambientValue, dt);

                for (int k=kstart; k<kend; k++)
                    for (int i=0; i<pConc.length; i++)
                        for (int j=0; j<pConc[0].length; j++) {
                            preSum += pConc[i][j][k];
                            boolean countCell = (mask == null) || mask[i][j][k];
                            if (countCell) {
                                preNot += pConc[i][j][k];
                            }
                            newValues[i][j][k] = pConc[i][j][k] + dt*(advect[i][j][k] + diffus[i][j][k]);
                            postSum += newValues[i][j][k];

                            /*if (newValues[i][j][k] < 0) {
                                System.out.println("negative result: i="+i+", j="+j+", k="+k);
                                System.out.println("counter = "+counter);
                                System.out.println("stop");
                            }*/

                            if (countCell) {
                                postNot += newValues[i][j][k];
                            }
                            if (feedingRate != null) {
                                if (k == FEEDING_KLEVEL)
                                    newValues[i][j][k] = newValues[i][j][k] + dt * feedingRateMult * feedingRate[i][j];
                            }

                            if (distFeeding != null) {
                                newValues[i][j][k] = newValues[i][j][k] + dt * feedingRateMult * distFeeding[i][j][k];
                            }

                        }
            }
        }

        int threads = Math.min(numProcessors,dim[2]);
        CalcPart[] cps = new CalcPart[threads];
        Thread[] thr = new Thread[threads];
        for (int i=0; i<threads; i++) {
            int kend = (dim[2]/threads)*(i+1);
            if (i == threads-1)
                kend = dim[2];

            cps[i] = new CalcPart((dim[2]/threads)*i, kend);
            thr[i] = new Thread(cps[i]);
            thr[i].start();
        }

        double preSum = 0, postSum = 0, preNot = 0, postNot = 0;

        try {
            for (int i=0; i<threads; i++) {
                thr[i].join();
                preSum += cps[i].preSum;
                postSum += cps[i].postSum;
                preNot += cps[i].preNot;
                postNot += cps[i].postNot;
            }

        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        // Copy new values into the state array:
        for (int i=0; i<newValues.length; i++)
            for (int j=0; j<newValues[i].length; j++)
            System.arraycopy(newValues[i][j], 0, pConc[i][j], 0, newValues[i][j].length);



        /*if (postSum > 1000 && postSum > preSum) {

            System.out.println("presum = "+preSum+", postsum = "+postSum);

            NetcdfFileWriteable ncfile = SaveNetCDF.initializeFile("c:\\work\\risingFeed\\presu_posu.nc", dim, 1, 1);
            SaveNetCDF.createCageVariables(ncfile, "presu", "posu");
            SaveNetCDF.saveCageVariable(ncfile, 0, "presu", presu, mask, true);
            SaveNetCDF.saveCageVariable(ncfile, 0, "posu", posu, mask, false);
            try {
                ncfile.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
            System.exit(0);
        }  */

        return new double[] {(-postSum+preSum)/dt, (-postNot+preNot)/dt};
    }

    /**
     * Activate linear reduction in ambient values dependent on distance measured along the normal
     * from the point ambRedCenter backwards along the vector dsVector.
     * @param varyAmbient
     */
    public void setVaryAmbient(boolean varyAmbient, /*double ambRedMult, double[] ambRedCenter, double[] dsVector,*/
                               double[] affinityProfile) {
        this.varyAmbient = varyAmbient;
        /*this.ambRedMult = ambRedMult;
        this.ambRedCenter = new double[2];
        this.ambRedCenter[0] = ambRedCenter[0];
        this.ambRedCenter[1] = ambRedCenter[1];
        this.dsVector = new double[2];
        this.dsVector[0] = dsVector[0];
        this.dsVector[1] = dsVector[1];*/
        this.affinityProfile = affinityProfile;

    }

    /**
     * Advection using a basic upstream scheme.
     */
    protected void calcAdvectAndDiff(int[] dim, int kstart, int kend, double[][][] feed, double dxy, double dz, boolean[][][] mask, double sinkingSpeed,
                                     double diffKappa, double[][][][] currentSpeed, double[] currentSpeedOffset, double[] ambientValue) {
        double[][][] nbr = new double[3][3][5];
        // Third dimension of nbr: above(diff), above(adv), same, below



        for (int k=kstart; k<kend; k++)
            for (int i=0; i<feed.length; i++)
                for (int j=0; j<feed[0].length; j++) {
                    // Skip if outside:
                    if ((mask != null) && !mask[i][j][k]) {
                        advect[i][j][k] = 0;
                        diffus[i][j][k] = 0;
                        continue;
                    }

                    // Get local currents. For each dimension we need the current on two edges:
                    double[][] current = new double[3][2];
                    current[0][0] = currentSpeed[i][j][k][0] + currentSpeedOffset[0];
                    current[1][0] = currentSpeed[i][j][k][1] + currentSpeedOffset[1];
                    current[2][0] = currentSpeed[i][j][k][2] + currentSpeedOffset[2] + sinkingSpeed;
                    current[0][1] = currentSpeed[i+1][j][k][0] + currentSpeedOffset[0];
                    current[1][1] = currentSpeed[i][j+1][k][1] + currentSpeedOffset[1];
                    current[2][1] = currentSpeed[i][j][k+1][2] + currentSpeedOffset[2] + sinkingSpeed;

                    // If feed is going upwards, stop at surface:
                    //if (k == 0 && current[2][0] < 0)
                    //    current[2][0] = 0;


                    //*** Defining cell neighbourhood:
                    nbr[1][1][2] = feed[i][j][k];

                    // Horizontal, -1 in x direction:
                    if (i == 0)
                        nbr[0][1][2] = ambientValue[k];
                    else {
                        nbr[0][1][2] = feed[i-1][j][k];
                        if ((mask != null) && !mask[i-1][j][k]) {
                            current[0][0] = 0;
                            nbr[0][1][2] = nbr[1][1][2];
                        }
                    }
                    // Horizontal, +1 in x direction:
                    //System.out.println("i = "+i);
                    //System.out.println("j = "+j);
                    //System.out.println("k = "+k);
                    if (i == dim[0]-1)
                        nbr[2][1][2] = ambientValue[k];
                    else {
                        nbr[2][1][2] = feed[i+1][j][k];
                        if ((mask != null) && !mask[i+1][j][k]) {
                            current[0][1] = 0;
                            nbr[2][1][2] = nbr[1][1][2];
                        }
                    }
                    // Horizontal, -1 in y direction:
                    if (j == 0)
                        nbr[1][0][2] = ambientValue[k];
                    else {
                        nbr[1][0][2] = feed[i][j-1][k];
                        if ((mask != null) && !mask[i][j-1][k]) {
                            current[1][0] = 0;
                            nbr[1][0][2] = nbr[1][1][2];
                        }
                    }
                    // Horizontal, +1 in y direction:
                    if (j == dim[1]-1)
                        nbr[1][2][2] = ambientValue[k];
                    else {
                        nbr[1][2][2] = feed[i][j+1][k];
                        if ((mask != null) && !mask[i][j+1][k]) {
                            current[1][1] = 0;
                            nbr[1][2][2] = nbr[1][1][2];
                        }
                    }

                    // Above:
                    if (k == 0) {
                        nbr[1][1][0] = feed[i][j][k]; // For diffusion
                        nbr[1][1][1] = 0.; // For advection
                    }
                    else {
                        nbr[1][1][0] = feed[i][j][k-1];
                        nbr[1][1][1] = feed[i][j][k-1];
                    }

                    // Below:
                    if (k == dim[2]-1) {
                        nbr[1][1][3] = ambientValue[k]; // By bottom
                        nbr[1][1][4] = ambientValue[k]; // By bottom
                    }
                    else {
                        nbr[1][1][3] = feed[i][j][k+1];
                        nbr[1][1][4] = feed[i][j][k+1];
                        /*if ((mask != null) && !mask[i][j][k+1]) {
                            current[2][1] = 0;
                            nbr[1][1][3] = nbr[1][1][2];
                        } */
                    }

                    // Note current directions:
                    int[][] curDir = new int[3][2];
                    curDir[0][0] = current[0][0] > 0 ? 1 : -1;
                    curDir[1][0] = current[1][0] > 0 ? 1 : -1;
                    curDir[2][0] = current[2][0] > 0 ? 1 : -1;
                    curDir[0][1] = current[0][1] > 0 ? 1 : -1;
                    curDir[1][1] = current[1][1] > 0 ? 1 : -1;
                    curDir[2][1] = current[2][1] > 0 ? 1 : -1;


                    /*if (k==25) {
                        System.out.println("current[2][0]="+current[2][0]+" current[2][1]"+current[2][1]);
                        System.out.println("nbr[1][1][1]="+nbr[1][1][1]+" nbr[1][1][2]="+nbr[1][1][2]+" nbr[1][1][4]="+nbr[1][1][4]);
                        System.out.println("ambientValue = "+ambientValue[k]);
                    } */
                    // Collect advection terms:
                    advect[i][j][k] = (1/dxy)*
                            (((current[0][0] > 0) ? current[0][0]*nbr[0][1][2] : current[0][0]*nbr[1][1][2])
                            +((current[0][1] < 0) ? -current[0][1]*nbr[2][1][2] : -current[0][1]*nbr[1][1][2])
                            +((current[1][0] > 0) ? current[1][0]*nbr[1][0][2] : current[1][0]*nbr[1][1][2])
                            +((current[1][1] < 0) ? -current[1][1]*nbr[1][2][2] : -current[1][1]*nbr[1][1][2]))
                           + (1/dz)*
                            (((current[2][0] > 0) ? current[2][0]*nbr[1][1][1] : current[2][0]*nbr[1][1][2])
                            +((current[2][1] < 0) ? -current[2][1]*nbr[1][1][4] : -current[2][1]*nbr[1][1][2]));





                    // Collect diffusion terms:
                    diffus[i][j][k] = diffKappa*((nbr[2][1][2] - 2*nbr[1][1][2] + nbr[0][1][2])/dxy/dxy
                        +(nbr[1][2][2] - 2*nbr[1][1][2] + nbr[1][0][2])/dxy/dxy
                        +(nbr[1][1][3] - 2*nbr[1][1][2] + nbr[1][1][0])/dz/dz);
                }

    }

    /** Advection using subgrid model with superbee slope limiter to reduce numerical diffusion:
     */
    protected void calcAdvectAndDiff2(int[] dim, int kstart, int kend, double[][][] feed, double dxy, double dz, boolean[][][] mask, double sinkingSpeed,
                                      double diffKappaXY, double diffKappaZ, double[][][][] currentSpeed, 
                                      double[] currentSpeedOffset, double[] ambientVals, double dt) {
        //double[][][] nbr = new double[5][5][7];
        boolean masking = mask != null;
        double c_h;
        double[] x_nb = new double[4], y_nb = new double[4], z_nb = new double[4],
                x_nb_diff = new double[2], y_nb_diff = new double[2], z_nb_diff = new double[2];

        double[] ambientValHere = new double[ambientVals.length];

        for (int k=kstart; k<kend; k++) {

            for (int i = 0; i < feed.length; i++)
                for (int j = 0; j < feed[0].length; j++) {
                    // Skip if outside:
                    if ((mask != null) && !mask[i][j][k]) {
                        advect[i][j][k] = 0;
                        diffus[i][j][k] = 0;
                        continue;
                    }

                    System.arraycopy(ambientVals, 0, ambientValHere, 0, ambientValHere.length);
                    if (varyAmbient) {
                        // Check if we are near an edge, otherwise the ambient value has no effect:
                        if (i<2 || j<2 || i>=feed.length-2 || j>=feed[0].length-2 || k<2 || k>=feed[0][0].length-2) {

                            //double[] vecHere = new double[]{((double) i) - ambRedCenter[0], ((double) j) - ambRedCenter[1]};
                            double ambRedValue = 0, ambRedValue2 = 0, multiplier = 0;
                            //double tmp1 = 0, tmp2 = 0;
                            double[] currHere = new double[] {currentSpeed[i][j][k][0] + currentSpeedOffset[0],
                                    currentSpeed[i][j][k][1] + currentSpeedOffset[1]};
                            double currHereSpeed = 100.*Math.sqrt(currHere[0]*currHere[0]+currHere[1]*currHere[1]);
                            //double multiplier = 1.1*Math.max(0.5, Math.min(2.5, 0.35*(7.4-currHereSpeed)));
                            double omega_red = 0.35;
                            if (currHereSpeed < 6)
                                omega_red = 1.85 - 0.25*currHereSpeed;

                            ambRedValue = ((double)(dim[0] - i)/((double)dim[0])) * omega_red * affinityProfile[k];



                            for (int ii = 0; ii < ambientValHere.length; ii++)
                                ambientValHere[ii] -= ambRedValue;

                        }


                    }
                    // Get local currents. For each dimension we need the current on two edges:
                    double[][] current = new double[3][2];
                    current[0][0] = currentSpeed[i][j][k][0] + currentSpeedOffset[0];
                    current[1][0] = currentSpeed[i][j][k][1] + currentSpeedOffset[1];
                    current[2][0] = currentSpeed[i][j][k][2] + sinkingSpeed + currentSpeedOffset[2];
                    current[0][1] = currentSpeed[i + 1][j][k][0] + currentSpeedOffset[0];
                    current[1][1] = currentSpeed[i][j + 1][k][1] + currentSpeedOffset[1];
                    current[2][1] = currentSpeed[i][j][k + 1][2] + sinkingSpeed + currentSpeedOffset[2];


                    //*** Defining cell neighbourhood:
                    c_h = feed[i][j][k];

                    // Horizontal, -1 in x direction:
                    x_nb[0] = pickVal(feed, ambientValHere, i - 2, j, k);
                    x_nb[1] = pickVal(feed, ambientValHere, i - 1, j, k);
                    // Check mask for diffusion:
                    if (!masking || checkMask(mask, false, i-1, j, k))
                        x_nb_diff[0] = x_nb[1];
                    else
                        x_nb_diff[0] = c_h;
                    // Check mask for advection:
                    if (i >= 0) {
                        if ((mask != null) && !mask[i - 1][j][k]) {
                            current[0][0] = 0;
                        }
                    }


                    // Horizontal, +1 in x direction:
                    x_nb[2] = pickVal(feed, ambientValHere, i + 1, j, k);
                    x_nb[3] = pickVal(feed, ambientValHere, i + 2, j, k);
                    // Check mask for diffusion:
                    if (!masking || checkMask(mask, false, i+1, j, k))
                        x_nb_diff[1] = x_nb[2];
                    else
                        x_nb_diff[1] = c_h;
                    // Check mask for advection:
                    if (i < dim[0] - 1) {
                        if ((mask != null) && !mask[i + 1][j][k]) {
                            current[0][1] = 0;
                        }
                    }
                    // Horizontal, -1 in y direction:
                    y_nb[0] = pickVal(feed, ambientValHere, i, j - 2, k);
                    y_nb[1] = pickVal(feed, ambientValHere, i, j - 1, k);
                    // Check mask for diffusion:
                    if (!masking || checkMask(mask, false, i, j-1, k))
                        y_nb_diff[0] = y_nb[1];
                    else
                        y_nb_diff[0] = c_h;
                    // Check mask for advection:
                    if (j > 0) {
                        if ((mask != null) && !mask[i][j - 1][k]) {
                            current[1][0] = 0;
                        }
                    }

                    // Horizontal, +1 in y direction:
                    y_nb[2] = pickVal(feed, ambientValHere, i, j + 1, k);
                    y_nb[3] = pickVal(feed, ambientValHere, i, j + 2, k);
                    // Check mask for diffusion:
                    if (!masking || checkMask(mask, false, i, j+1, k))
                        y_nb_diff[1] = y_nb[2];
                    else
                        y_nb_diff[1] = c_h;
                    // Check mask for advection:
                    if (j < dim[1] - 1) {
                        if ((mask != null) && !mask[i][j + 1][k]) {
                            current[1][1] = 0;
                        }
                    }

                    // Above:
                    z_nb[0] = pickVal(feed, ambientValHere, i, j, k - 2);
                    z_nb[1] = pickVal(feed, ambientValHere, i, j, k - 1);
                    z_nb_diff[0] = z_nb[1];
                    if (k == 0)
                        z_nb_diff[0] = c_h; // For diffusion

                    // Below:
                    z_nb[2] = pickVal(feed, ambientValHere, i, j, k + 1);
                    z_nb[3] = pickVal(feed, ambientValHere, i, j, k + 2);
                    z_nb_diff[1] = z_nb[2];
                    // Check mask for diffusion:
                    if (masking && !checkMask(mask, false, i, j, k+1))
                        z_nb_diff[0] = c_h;

                    // Note current directions:
                    int[][] curDir = new int[3][2];
                    curDir[0][0] = current[0][0] > 0 ? 1 : -1;
                    curDir[1][0] = current[1][0] > 0 ? 1 : -1;
                    curDir[2][0] = current[2][0] > 0 ? 1 : -1;
                    curDir[0][1] = current[0][1] > 0 ? 1 : -1;
                    curDir[1][1] = current[1][1] > 0 ? 1 : -1;
                    curDir[2][1] = current[2][1] > 0 ? 1 : -1;


                    boolean simpleAdv = false;

                    if (!simpleAdv) {
                        // Collect advection terms:
                        advect[i][j][k] = superbeeAdv(dt, dxy, x_nb[0], x_nb[1], c_h, x_nb[2], x_nb[3], current[0][0], current[0][1])
                                + superbeeAdv(dt, dxy, y_nb[0], y_nb[1], c_h, y_nb[2], y_nb[3], current[1][0], current[1][1])
                                + superbeeAdv(dt, dz, z_nb[0], z_nb[1], c_h, z_nb[2], z_nb[3], current[2][0], current[2][1]);
                    }
                    else {

                        advect[i][j][k] = simpleAdv(dt, dxy, x_nb[0], x_nb[1], c_h, x_nb[2], x_nb[3], current[0][0], current[0][1])
                                + simpleAdv(dt, dxy, y_nb[0], y_nb[1], c_h, y_nb[2], y_nb[3], current[1][0], current[1][1])
                                + simpleAdv(dt, dz, z_nb[0], z_nb[1], c_h, z_nb[2], z_nb[3], current[2][0], current[2][1]);

                        /*advect[i][j][k] = simpleAdv2(dt, dxy, x_nb[0], x_nb[1], c_h, x_nb[2], x_nb[3], current[0][0], current[0][1])
                                + simpleAdv2(dt, dxy, y_nb[0], y_nb[1], c_h, y_nb[2], y_nb[3], current[1][0], current[1][1])
                                + simpleAdv2(dt, dz, z_nb[0], z_nb[1], c_h, z_nb[2], z_nb[3], current[2][0], current[2][1]);*/

                    }

                    // Collect diffusion terms:
                    diffus[i][j][k] = diffKappaXY * ((x_nb_diff[0] - 2 * c_h + x_nb_diff[1]) / dxy / dxy
                            + (y_nb_diff[0] - 2 * c_h + y_nb_diff[1]) / dxy / dxy)
                            + diffKappaZ * ((z_nb_diff[0] - 2 * c_h + z_nb_diff[1]) / dz / dz);


                }
        }
    }


    private double pickVal(double[][][] fc, double[] ambVal, int i, int j, int k) {
        if (k < 0) return 0;
        if (i<0 || i>=fc.length || j<0 || j>=fc[0].length || k>=fc[0][0].length)
            return ambVal[Math.min(ambVal.length-1,k)];
        else return fc[i][j][k];
    }

    /**
     * Return the value of the mask at the given coordinate, or the default value if coordinates are outside the domain.
     */
    private boolean checkMask(boolean[][][] mask, boolean defVal, int i, int j, int k) {
        if (i<0 || i>=mask.length || j<0 || j>=mask[0].length || k<0 || k>=mask[0][0].length)
            return defVal;
        else return mask[i][j][k];
    }

    private double minModAdv(double dt, double dx, double c_ll, double c_l, double c_c, double c_r, double c_rr, double v_l, double v_r) {
        double sum = 0;
        if (v_l >= 0) {
            double sigma_l = minmod((c_c-c_l)/dx, (c_l-c_ll)/dx);
            sum += (v_l/dx)*(c_l + (sigma_l/2.0)*(dx-v_l*dt));
        }
        else {
            double sigma_l = minmod((c_c-c_l)/dx, (c_r-c_c)/dx);
            sum += (v_l/dx)*(c_c - (sigma_l/2.0)*(dx+v_l*dt));
        }

        if (v_r >= 0) {
            double sigma_r = minmod((c_c-c_l)/dx, (c_r-c_c)/dx);
            sum -= (v_r/dx)*(c_c + (sigma_r/2.0)*(dx-v_r*dt));
        }
        else {
            double sigma_r = minmod((c_r-c_c)/dx, (c_rr-c_r)/dx);
            sum -= (v_r/dx)*(c_r - (sigma_r/2.0)*(dx+v_r*dt));
        }

        return sum;
    }


    private double simpleAdv(double dt, double dx, double c_ll, double c_l, double c_c, double c_r, double c_rr, double v_l, double v_r) {
        double sum = 0;
        if (v_l > 0) {
            sum += (c_l - c_c) * v_l;
        }
        if (v_r < 0) {
            sum -= (c_r - c_c) * v_r;
        }
        return sum;
    }

    private double simpleAdv2(double dt, double dx, double c_ll, double c_l, double c_c, double c_r, double c_rr, double v_l, double v_r) {
        double sum = 0;
        if (v_l > 0) {
            sum += c_l * v_l;
        } else {
            sum += c_c * v_l;
        }
        if (v_r < 0) {
            sum -= c_r * v_r;
        } else {
            sum -= c_c * v_r;
        }
        return sum;
    }

    private double superbeeAdv(double dt, double dx, double c_ll, double c_l, double c_c, double c_r, double c_rr, double v_l, double v_r) {
        double sum = 0;
        if (v_l >= 0) {
            double sigma_l = maxmod(minmod((c_c-c_l)/dx, 2.*(c_l-c_ll)/dx),
                        minmod(2.*(c_c-c_l)/dx, (c_l-c_ll)/dx));
            sum += (v_l/dx)*(c_l + (sigma_l/2.0)*(dx-v_l*dt));
        }
        else {
            double sigma_l = maxmod(minmod((c_c-c_l)/dx, 2.*(c_r-c_c)/dx),
                        minmod(2.*(c_c-c_l)/dx, (c_r-c_c)/dx));
            sum += (v_l/dx)*(c_c - (sigma_l/2.0)*(dx+v_l*dt));
        }

        if (v_r >= 0) {
            double sigma_r = maxmod(minmod((c_c-c_l)/dx, 2.*(c_r-c_c)/dx),
                        minmod(2.*(c_c-c_l)/dx, (c_r-c_c)/dx));
            sum -= (v_r/dx)*(c_c + (sigma_r/2.0)*(dx-v_r*dt));

        }
        else {
            double sigma_r = maxmod(minmod((c_r-c_c)/dx, 2.*(c_rr-c_r)/dx),
                        minmod(2.*(c_r-c_c)/dx, (c_rr-c_r)/dx));
            sum -= (v_r/dx)*(c_r - (sigma_r/2.0)*(dx+v_r*dt));
        }

        return sum;
    }

    private double maxmod(double a, double b) {
        if (a*b < 0)
            return 0;
        else return (Math.abs(a) > Math.abs(b)) ? a : b;
    }

    private double minmod(double a, double b) {
        if (a*b < 0)
            return 0;
        else return (Math.abs(a) < Math.abs(b)) ? a : b;
    }


    public static double getMeanFeedDepth(double[][][] feed, double dz, boolean[][][] mask) {
        double totFeed = 0, weightedSum = 0;
        for (int i=0; i<feed.length; i++)
            for (int j=0; j<feed[0].length; j++)
                for (int k=0; k<feed[0][0].length; k++) {
                    if ((mask != null) && !mask[i][j][k])
                        continue;
                    totFeed += feed[i][j][k];
                    weightedSum += feed[i][j][k]*((double)k + 0.5)*dz;
                }
        if (totFeed > 0)
            return weightedSum/totFeed;
        else return Double.NaN;
    }

    public void clearOutside(double[][][] feed, boolean[][][] mask) {
        for (int i=0; i<feed.length; i++)
            for (int j=0; j<feed[0].length; j++)
                for (int k=0; k<feed[0][0].length; k++) {
                    if ((mask != null) && !mask[i][j][k])
                        feed[i][j][k] = 0;
                }

    }

    /**
     * Create a Gaussian distribution of feeding input around a single cell with a given standard deviation.
     * The distribution is normalized so the sum of input values is 1.
     */
    public static void setNormalFeedingDist(double[][] feedingDist, int centerX, int centerY, double dxy, double sigma) {
        double sum = 0;
        for (int i=0; i<feedingDist.length; i++)
            for (int j=0; j<feedingDist[0].length; j++) {
                double sqDistance = Math.pow((double)(i - centerX)*dxy, 2.0) + Math.pow((double)(j - centerY)*dxy, 2.0);
                feedingDist[i][j] = Math.exp(-sqDistance/(2*sigma*sigma));
                sum += feedingDist[i][j];
            }
        double sum2 = 0;
        for (int i=0; i<feedingDist.length; i++)
            for (int j=0; j<feedingDist[0].length; j++) {
                feedingDist[i][j] /= sum;
                sum2 += feedingDist[i][j];
            }

    }

    public static void setUniformFeedingDist(double[][] feedingDist, int padding) {
        double sum = 0;
        for (int i=0; i<feedingDist.length; i++)
            for (int j=0; j<feedingDist[0].length; j++) {
                if ((i >= padding) && (i < feedingDist.length-padding) && (j >= padding) && (j < feedingDist[0].length-padding))
                    feedingDist[i][j] = 1;
                else
                    feedingDist[i][j] = 0;
                sum += feedingDist[i][j];
            }
        double sum2 = 0;
        for (int i=0; i<feedingDist.length; i++)
            for (int j=0; j<feedingDist[0].length; j++) {
                feedingDist[i][j] /= sum;
                sum2 += feedingDist[i][j];
            }

    }



    public static void setSinglePointFeedingDist(double[][] feedingDist, int[] pos) {
        double sum = 0;
        for (int i=0; i<feedingDist.length; i++)
            for (int j=0; j<feedingDist[0].length; j++) {
                if ((i == pos[0]) && (j == pos[1]))
                    feedingDist[i][j] = 1;
                else
                    feedingDist[i][j] = 0;
                sum += feedingDist[i][j];
            }
        double sum2 = 0;
        for (int i=0; i<feedingDist.length; i++)
            for (int j=0; j<feedingDist[0].length; j++) {
                feedingDist[i][j] /= sum;
                sum2 += feedingDist[i][j];
            }

    }

    public static void initField(double[][][] field, double value) {
        for (int i=0; i<field.length; i++)
            for (int j=0; j<field[i].length; j++)
                for (int k=0; k<field[i][j].length; k++)
                    field[i][j][k] = value;
    }

}