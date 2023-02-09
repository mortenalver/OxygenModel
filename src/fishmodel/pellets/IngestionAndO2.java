package fishmodel.pellets;

/**
 * Calculate pellet ingestion by fish according to the model of Alver et al. (2004).
 */
public class IngestionAndO2 {

    final static double
        T_h = 12, // Handling time
        k_T_s = 1., // Factor for what feed particle count related to fish count makes search time important
        b = 0.4, // Exponent for confusion factor
        c = 0.5; // Exponent for f_d factor


    final static double o2consumptionMult = 1.3; // This factor can be used to globally multiply the o2 consumption of the fish.

    static double U = 1; // Swimming speed (body lengths/s)

    // We calculate O2 consumption according to a hybrid pattern consisting of two parts:
    // 1. Proportionally to feed distribution. This component is only present when fish is feeding
    // 2. Even distribution over cage volume.
    // When fish is feeding, the factor o2_even_fraction determines the fraction that is determined
    // according to part 2.
    static double o2_even_fraction = 0.75;//0.5;

    // O2 consumption (Grøttum and Sigholt, 1998):
    //    VO2 (mg/kg/h) = 61.6 * BW^0.33 * 1.03^T * 1.79^U
    // BW: body weight (kg)
    // T: temperature  (C)
    // U: swimming speed (body lengths/s)

    public static double[] calculateIngestion(double dt, double[][][] feed, double[][][] o2, double[][][] affinity, double[][][] o2Affinity, double o2AffSum,
                                              double[][][] ingDist, double dxy, double dz, boolean[][][] mask, double pelletWeight, double T_w, SimpleFish fish) {
        double N = fish.getTotalN();
        double WtotKg = 0.001*fish.getTotalW();
        if (N == 0)
            return new double[] {0, 0};


        double totalFeed = 0;
        for (int i=0; i<feed.length; i++)
            for (int j=0; j<feed[0].length; j++)
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
            for (int i=0; i<feed.length; i++)
                for (int j=0; j<feed[0].length; j++)
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

        // TODO: need to find good way to model effect of o2 level on appetite
        //double o2Thresh = -70;
        //double o2AppetiteMult = 1.0;
        /*if (o2PercAverage < o2Thresh)
            o2AppetiteMult = Math.max(0, o2PercAverage/o2Thresh);*/

        //System.out.println("o2PercAverage = "+o2PercAverage);
        //System.out.println("MaxDensity = "+maxDensity);

        /*
        // OLD rho calculation:
        double c_ss = 0.;
        double nCells = 0;
        if (feeding)
            for (int i=0; i<feed.length; i++)
                for (int j=0; j<feed[0].length; j++)
                    for (int k=0; k<feed[0][0].length; k++)
                        if ((mask == null) || mask[i][j][k]) {
                            c_ss += affinity[i][j][k]*Math.pow(feed[i][j][k]/totalFeed, 2.0);
                            nCells += affinity[i][j][k];
                        }
        rho = feeding ? 1/(nCells*c_ss) : 0;*/

        // Confusion factor:
        double p_c = Math.pow(rho, b);

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
        //System.out.println("Total intake = "+totalIntake);
        if (totalIntake == 0)
            feeding = false;
        //System.out.println("Relative removal = "+totalIntake*dt/totalFeed);
        // Make sure the ingestion doesn't exceed the available feed:
        if (dt*totalIntake > totalFeed) {
            double multiplier = totalFeed/(dt*totalIntake);
            totalIntake *= multiplier;
            for (int i=0; i<fish.getNGroups(); i++)
                w_f[i] *= multiplier;
        }

        // Calculate total oxygen consumption in mg per second:
        double o2Cons = 0;
        if (o2 != null) {
            for (int i=0; i<fish.getNGroups(); i++) {
                // Consumption: N multiplied by weight in kg multiplied by O2 consumption per kg divided by 3600 s:
                o2Cons += o2consumptionMult*fish.getN(i)*0.001*fish.getW(i)*61.6*Math.pow(fish.getW(i)*0.001, -0.33)*Math.pow(1.03, T_w)*Math.pow(1.79, U)/3600.0;
                //System.out.println(i+": \t"+(/*0.001*fish.getW(i)**/61.6*Math.pow(fish.getW(i)*0.001, -0.33)*Math.pow(1.03, T_w)*Math.pow(1.79, U)));
            }
            //System.out.println("o2: "+o2Cons);
        }

        // Calculate how much is removed from each cell:
        double[][][] cellIng = new double[feed.length][feed[0].length][feed[0][0].length];
        double sumCellIng = 0;
        double correction = 1;
        if (feeding) {
            for (int i=0; i<feed.length; i++)
                for (int j=0; j<feed[0].length; j++)
                    for (int k=0; k<feed[0][0].length; k++) {
                        // Remove same relative fraction of feed everywhere:
                        cellIng[i][j][k] = affinity[i][j][k]*totalIntake*feed[i][j][k]/totalFeed;
                        sumCellIng += cellIng[i][j][k];


                    }
            correction = sumCellIng > 0 ? totalIntake/sumCellIng : 0;
        }


        // Find the O2 consumption in the two parts:
        double o2Cons_even = feeding ? o2_even_fraction*o2Cons : o2Cons;
        double o2Cons_feeding = feeding ? (1.-o2_even_fraction)*o2Cons : 0.;

        // Remove feed and o2 from cells:
        double presum = 0, postsum = 0;
        for (int i=0; i<feed.length; i++)
            for (int j=0; j<feed[0].length; j++)
                for (int k=0; k<feed[0][0].length; k++) {
                    if ((mask == null) || mask[i][j][k]) {

                        feed[i][j][k] -= dt*correction*cellIng[i][j][k];
                        if (ingDist != null) ingDist[i][j][k] = correction*cellIng[i][j][k];

                        if (o2 != null) {
                            presum += o2[i][j][k];
                            // TODO: negative o2 values are simply cut off, no reduction of consumption when o2 is low
                            // Since O2 is given as a concentration (mg/l), we need to divide by the cell volume in l:
                            // Remove part of the O2 consumption according to the same distribution as feed ingestion.
                            if (feeding)
                                o2[i][j][k] = Math.max(0., o2[i][j][k] - dt*(correction*cellIng[i][j][k]/totalIntake)*o2Cons_feeding/(1000.0*dxy*dxy*dz));

                            // Remove part (or all) of O2 consumption according to affinity:
                            o2[i][j][k] = Math.max(0., o2[i][j][k] - dt*(o2Affinity[i][j][k]/o2AffSum)*o2Cons_even/(1000.0*dxy*dxy*dz));

                            postsum += o2[i][j][k];

                        }
                    }
                }

        double o2ConsumptionRate = (presum-postsum)*1000.0*dxy*dxy*dz/dt; // mg o2 removed from volume per second

        // Add feed to stomachs:
        for (int i=0; i<fish.getNGroups(); i++) {
            fish.stepGutContent(i, dt, T_w);
            fish.addIngestion(i, dt*w_f[i]);
            fish.setIngRate(i, w_f[i]);
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
