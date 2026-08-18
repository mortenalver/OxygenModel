package fishmodel.pellets;


/**
 * Calculate pellet ingestion by fish according to the model of Alver et al. (2004).
 * This variant operates on a subgrid representing one cage in a multi-cage setup
 */
public class IngestionAndO2Subgrid {


    final static double
        T_h = 12, // Handling time
        k_T_s = 1., // Factor for what feed particle count related to fish count makes search time important
        b = 0.4, // Exponent for confusion factor
        c = 0.5; // Exponent for f_d factor

    // This factor can be used to globally multiply the o2 consumption of the fish when using the Grøttum & Sigholt model:
    public static double o2consumptionMultOld = 1.2*1.3;

    // This factor can be used to globally multiply the o2 consumption of the fish when using the new o2 consumption model:
    public static double o2consumptionMultNew = 1.0;

    private static boolean addDigestiveO2Cons = false;

    static double U = 1; // Swimming speed (body lengths/s)

    // We calculate O2 consumption according to a hybrid pattern consisting of two parts:
    // 1. Proportionally to feed distribution. This component is only present when fish is feeding
    // 2. Even distribution over cage volume.
    // When fish is feeding, the factor o2_even_fraction determines the fraction that is determined
    // according to part 2.
    // Ocean Farm 1: 0.75
    static double o2_even_fraction = 0.75;//0.5;

    public static void setAddDigestiveO2Cons(boolean value) {
        IngestionAndO2Subgrid.addDigestiveO2Cons = value;
    }

    public static void setO2EvenFraction(double value) { IngestionAndO2Subgrid.o2_even_fraction = value; }

    // O2 consumption (Grøttum and Sigholt, 1998):
    //    VO2 (mg/kg/h) = 61.6 * BW^0.33 * 1.03^T * 1.79^U
    // BW: body weight (kg)
    // T: temperature  (C)
    // U: swimming speed (body lengths/s)

    public static double[] calculateIngestion(double dt, double[][][] feed, double[][][] o2, double[][][] affinity, double[][][] o2Affinity, double o2AffSum,
                                              int[][] ranges,
                                              double[][][] ingDist, double[][][] o2ConsDist, double dxy, double dz, boolean[][][] mask, double pelletWeight,
                                              double[] T_w, SimpleFish fish, double o2Cons_perturb, boolean useNewConsumptionModel) {
        double N = fish.getTotalN();
        double WtotKg = 0.001*fish.getTotalW();
        //System.out.println("N="+N+" / totW="+WtotKg);
        if (N == 0)
            return new double[] {0, 0, 0};

        // x and y ranges for calculation:
        int x0 = ranges[0][0],
                x1 = ranges[0][1],
                y0 = ranges[1][0],
                y1 = ranges[1][1];

        double totalFeed = 0;
        for (int i=x0; i<x1; i++)
            for (int j=y0; j<y1; j++)
                for (int k=0; k<feed[0][0].length; k++)
                    if ((mask == null) || mask[i][j][k]) {
                        totalFeed += affinity[i][j][k]*feed[i][j][k];
                    }

        boolean feeding = totalFeed > 1e-3;

        double w_0 = feeding ? 1./(T_h + k_T_s*N*pelletWeight/totalFeed) : 0;

        // NEW rho calculation:
        double muSum = 0;
        double wsum = 0;
        double cellVol = dxy*dxy*dz; // Cell volume in m3
        double maxDensity = 0;
        if (feeding)
            for (int i=x0; i<x1; i++)
                for (int j=y0; j<y1; j++)
                    for (int k=0; k<feed[0][0].length; k++)
                        if ((mask == null) || mask[i][j][k]) {
                            double wHere = WtotKg*affinity[i][j][k]*feed[i][j][k]/totalFeed;
                            wsum += wHere;
                            double density = wHere/cellVol;
                            if (density > maxDensity)
                                maxDensity = density;
                            muSum += wHere*mu(density);


                        }
        double rho = feeding ? muSum/WtotKg : 0;

        // Confusion factor:
        double p_c = rho > 0 ? Math.pow(rho, b) : 0;


        double[] p_a = new double[fish.getNGroups()];
        double totalW = fish.getTotalW();
        double f_a = 0;
        // Calculate appetite factors and the f_a factor:
        for (int i=0; i<fish.getNGroups(); i++) {
            p_a[i] = fish.getAppetite(i);
            f_a += fish.getN(i)*p_a[i]*fish.getW(i)/totalW;
        }
        //System.out.println("f_a = "+f_a);
        double f_d = rho > 0 ? Math.pow(rho, -c) : 0;

        // Calculate p_h factors and feed intake per group:
        double maxW = fish.getMaxW();
        double[] w_f = new double[fish.getNGroups()];
        double totalIntake = 0;
        for (int i=0; i<fish.getNGroups(); i++) {

            // Hierarchy factor:
            double p_h = Math.pow(fish.getW(i)/maxW, f_a*f_d);

            w_f[i] = pelletWeight*w_0*p_c*p_a[i]*p_h;
            totalIntake += fish.getN(i)*w_f[i];

            //System.out.println("w_f["+i+"] = "+w_f[i]);
        }

        if (totalIntake == 0)
            feeding = false;
        //System.out.println("Relative removal = "+totalIntake*dt/totalFeed);
        // Make sure the ingestion doesn't exceed the available feed:
        if (dt*totalIntake > totalFeed) {
            double multiplier = totalIntake > 0 ? totalFeed / (dt * totalIntake) : 0;
            totalIntake *= multiplier;
            //System.out.println("Multiplier: "+multiplier);
            for (int i = 0; i < fish.getNGroups(); i++) {
                w_f[i] *= multiplier;

            }
        }

        // Calculate how much is removed from each cell:
        double[][][] cellIng = new double[feed.length][feed[0].length][feed[0][0].length];
        double sumCellIng = 0;
        double correction = 1;
        double finalScaling = 1.0;
        if (feeding) {
            for (int i = x0; i < x1; i++)
                for (int j = y0; j < y1; j++)
                    for (int k = 0; k < feed[0][0].length; k++) {
                        // Remove same relative fraction of feed everywhere:
                        cellIng[i][j][k] = affinity[i][j][k] * totalIntake * feed[i][j][k] / totalFeed;
                        sumCellIng += cellIng[i][j][k];
                        //System.out.println(affinity[i][j][k]);
                    }
            correction = sumCellIng > 0 ? totalIntake / sumCellIng : 0;
            //System.out.println("Correction: "+correction);
            // Remove feed from cells:
            double removed = 0;
            for (int i = x0; i < x1; i++)
                for (int j = y0; j < y1; j++)
                    for (int k = 0; k < feed[0][0].length; k++) {
                        if ((mask == null) || mask[i][j][k]) {
                            double toRemove = Math.min(feed[i][j][k], dt * correction * cellIng[i][j][k]);
                            feed[i][j][k] -= toRemove;
                            removed += toRemove;
                            if (ingDist != null) ingDist[i][j][k] = correction * cellIng[i][j][k];
                        }

                        if (Double.isNaN(feed[i][j][k])) {
                            System.out.println("NaN");
                        }
                    }

            //System.out.println("Relative feed removal success: "+removed/(dt*correction*totalIntake));
            if (removed < dt*correction*totalIntake) {
                finalScaling = removed/(dt*correction*totalIntake);
                totalIntake *= finalScaling;
            }
        }


        // Oxygen consumption.
        // Step 1: compute the affinity-dependent distribution of the fish with regard to oxygen.
        // betaBar should sum up to 1.0, however this doesn't hold when looking at a subgrid. Therefore
        // we use 1/sumbb as a correction factor when we use the betaBar values.
        // Each element gives the fraction of O2 ingesting fish in one cell
        double sumbb = 0; // Sum of all betaBar values. Its inverse is used as correction factor.
        double[][][] betaBar = new double[o2.length][o2[0].length][o2[0][0].length];
        for (int i=x0; i<x1; i++)
            for (int j=y0; j<y1; j++)
                for (int k=0; k<feed[0][0].length; k++) {

                    if ((mask == null) || mask[i][j][k]) {
                        if (feeding) {
                            betaBar[i][j][k] = (correction * cellIng[i][j][k] / totalIntake) * (1. - o2_even_fraction);

                            // Remove part (or all) of O2 consumption according to affinity:
                            betaBar[i][j][k] += (o2Affinity[i][j][k] / o2AffSum) * o2_even_fraction;
                        }
                        else {
                            // Remove part (or all) of O2 consumption according to affinity:
                            betaBar[i][j][k] += (o2Affinity[i][j][k] / o2AffSum);
                        }
                    }
                    sumbb += betaBar[i][j][k];
                }
        //System.out.println("sum betabar = "+sumbb);

        // Step 2: cycle through all cells, and compute the O2 consumption given the amount of fish and temperature
        // in that cell
        double presum=0., postsum = 0.;
        for (int i=x0; i<x1; i++)
            for (int j=y0; j<y1; j++)
                for (int k=0; k<feed[0][0].length; k++) {
                    if ((mask == null) || mask[i][j][k]) {
                        double consHere = 0;
                        for (int kg=0; kg< fish.getNGroups(); kg++) {
                            //consHere += o2consumptionMult*(1.0 + o2Cons_perturb)*fish.getN(kg)*0.001*fish.getW(kg)*61.6*Math.pow(fish.getW(kg)
                            //        *0.001, -0.33)*Math.pow(1.03, T_w[k])*Math.pow(1.79, U)/3600.0;

                            double o2ConsumptionGandS = 0;

                            if (!useNewConsumptionModel) {
                                // Grøttum and Sigholt:
                                o2ConsumptionGandS = o2consumptionMultOld*fish.getN(kg)*0.001*fish.getW(kg)*61.6*Math.pow(fish.getW(kg)
                                        *0.001, -0.33)*Math.pow(1.03, T_w[k])*Math.pow(1.79, U)/3600.0;
                            } else {

                                // Orig submitted: 93.9*W^-0.13*1.03.^T*1.56.^U
                                // Final version: 79.7*W^-0.14*1.04^T*1.64^U
                                o2ConsumptionGandS = o2consumptionMultNew*fish.getN(kg)*0.001*fish.getW(kg)*79.7*Math.pow(fish.getW(kg)
                                        *0.001, -0.14)*Math.pow(1.04, T_w[k])*Math.pow(1.64, U)/3600.0;
                            }

                            if (addDigestiveO2Cons) {
                                // Get rate of O2 consumption from digestion in g/s per individual. Multiply by N:
                                double digO2Cons = fish.getN(kg)*fish.calcDigestiveO2Cons(kg, T_w[k]);

                                o2ConsumptionGandS += digO2Cons;
                            }

                            consHere += (1.0 + o2Cons_perturb)*o2ConsumptionGandS;
                        }
                        consHere *= betaBar[i][j][k]/sumbb;
                        presum += o2[i][j][k];
                        // TODO: negative o2 values are simply cut off, no reduction of consumption when o2 is low
                        // Since O2 is given as a concentration (mg/l), we need to divide by the cell volume in l:
                        o2[i][j][k] = Math.max(0., o2[i][j][k] - dt*consHere/(1000.0*dxy*dxy*dz));
                        if (o2ConsDist != null) {
                            o2ConsDist[i][j][k] = consHere/(1000.0*dxy*dxy*dz);
                        }
                        postsum += o2[i][j][k];
                    }
                }


        double o2ConsumptionRate = (presum-postsum)*1000.0*dxy*dxy*dz/dt; // mg o2 removed from volume per second

        //System.out.println("consSum="+consSum+" , consRate="+o2ConsumptionRate);

        // Add feed to stomachs:
        for (int i=0; i<fish.getNGroups(); i++) {
            fish.stepGutContent(i, dt, T_w[0]);
            fish.addIngestion(i, dt*w_f[i]*finalScaling);
            fish.setIngRate(i, w_f[i]*finalScaling);
        }

        return new double[] {totalIntake, rho, o2ConsumptionRate};
    }

    /**
     * Calculate local mu value for a given fish density (kg/m3)
     * @param density
     * @return
     */
    public static double mu(double density) {
        double thresh = 50;//110;
        if (density < thresh) return 1;
        else return Math.max(0, 1-(density-thresh)/50);

    }
}
