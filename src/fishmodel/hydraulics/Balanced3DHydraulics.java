package fishmodel.hydraulics;

import java.text.NumberFormat;
import java.util.Locale;
import java.util.Random;

public class Balanced3DHydraulics {

    public static void getTurbulentHydraulicField(int[] dims, double dx, double rad, double currRedFactor,
                                                double[][] profile, double[][][][] uvw) {
        // Preliminaries:
        Random rand = new Random();
        double centerPosX = 0.5*(dims[0]-1),
                centerPosY = 0.5*(dims[1]-1),
                turbStd = 0.03;
        int gaussianBlurRadius = 2;

        // Set up u values:
        double[][][] uu = new double[dims[0]+1][dims[1]+1][dims[2]];
        for (int i=0; i<dims[0]+1; i++)
            for (int j=0; j<dims[1]; j++)
                for (int k=0; k<dims[2]; k++) {
                    // Cell center horizontal distance from domain center:
                    double cvecX = ((double)i)-centerPosX,
                            cvecY = ((double)j + 0.5)-centerPosY;
                    double normX = cvecY, normY = -cvecX; // Normal vector to the one from the center
                    double cdist = Math.sqrt(cvecX*cvecX + cvecY*cvecY);
                    uu[i][j][k] = profile[k][0]*currRedFactor + rand.nextGaussian()*turbStd;


                }
        // Smooth array:
        uu = ArraySmoothing3D.applyGaussianBlur3D(uu, 2);
        // Copy to output array:
        for (int i=0; i<dims[0]+1; i++)
            for (int j=0; j<dims[1]; j++)
                for (int k=0; k<dims[2]; k++) {
                    uvw[i][j][k][0] = uu[i][j][k];
                }

        // Set up v values:
        double[][][] vv = new double[dims[0]+1][dims[1]+1][dims[2]];
        for (int i=0; i<dims[0]; i++)
            for (int j=0; j<dims[1]+1; j++)
                for (int k=0; k<dims[2]; k++) {
                    // Cell center horizontal distance from domain center:
                    double cvecX = ((double)i)-centerPosX,
                            cvecY = ((double)j + 0.5)-centerPosY;
                    double normX = cvecY, normY = -cvecX; // Normal vector to the one from the center
                    double cdist = Math.sqrt(cvecX*cvecX + cvecY*cvecY);
                    vv[i][j][k] = profile[k][1]*currRedFactor + rand.nextGaussian()*turbStd;


                }
        // Smooth array:
        vv = ArraySmoothing3D.applyGaussianBlur3D(vv, 2);
        // Copy to output array:
        for (int i=0; i<dims[0]; i++)
            for (int j=0; j<dims[1]+1; j++)
                for (int k=0; k<dims[2]; k++) {
                    uvw[i][j][k][1] = vv[i][j][k];
                }

        // Set up w values to balance the flow field in all cells:
        computeVerticalSpeeds(dims, uvw);
    }

    public static void getBjoroyaHydraulicField(int[] dims, double dx, double rad, double currRedFactor,
                                                double[][] profile, double[][][][] uvw) {
        // Preliminaries:
        boolean outFlow = false;
        boolean circularFlow = true;
        double centerPosX = 0.5*(dims[0]-1),
                centerPosY = 0.5*(dims[1]-1),
                gamma = 0.5, // Sharpness of current reduction at cage wall
                outFlowRate = 0*0.01,
                outFlowMaxLayer = 2,
                circularFlowRate = -0.01;

        // Set up u values:
        for (int i=0; i<dims[0]+1; i++)
            for (int j=0; j<dims[1]; j++)
                for (int k=0; k<dims[2]; k++) {
                    // Cell center horizontal distance from domain center:
                    double cvecX = ((double)i)-centerPosX,
                            cvecY = ((double)j + 0.5)-centerPosY;
                    double normX = cvecY, normY = -cvecX; // Normal vector to the one from the center
                    double cdist = Math.sqrt(cvecX*cvecX + cvecY*cvecY);
                    uvw[i][j][k][0] = profile[k][0]*(1-currRedFactor*sigmoidFun(rad, cdist*dx, gamma));

                    if (outFlow && k<outFlowMaxLayer && cdist > 0) {
                        uvw[i][j][k][0] += outFlowRate * sigmoidFun(rad, cdist*dx, gamma) * cvecX / cdist;
                    }

                    if (circularFlow && cdist*dx <= rad) {
                        double cfr = circularFlowRate * (cdist*dx/rad - Math.pow(cdist*dx/rad, 3.0));
                        uvw[i][j][k][0] += normX*cfr;
                    }
                }

        // Set up v values:
        for (int i=0; i<dims[0]; i++)
            for (int j=0; j<dims[1]+1; j++)
                for (int k=0; k<dims[2]; k++) {
                    // Cell center horizontal distance from domain center:
                    double cvecX = ((double)i + 0.5)-centerPosX,
                            cvecY = ((double)j)-centerPosY;
                    double normX = cvecY, normY = -cvecX; // Normal vector to the one from the center
                    double cdist = Math.sqrt(cvecX*cvecX + cvecY*cvecY);
                    uvw[i][j][k][1] = profile[k][1]*(1-currRedFactor*sigmoidFun(rad, cdist*dx, gamma));

                    if (outFlow && k<outFlowMaxLayer && cdist > 0) {
                        uvw[i][j][k][1] += outFlowRate * sigmoidFun(rad, cdist*dx, gamma) * cvecY / cdist;
                    }

                    if (circularFlow && cdist*dx <= rad) {
                        double cfr = circularFlowRate * (cdist*dx/rad - Math.pow(cdist*dx/rad, 3.0));
                        uvw[i][j][k][1] += normY*cfr;
                    }
                }

        // Set up w values to balance the flow field in all cells:
        computeVerticalSpeeds(dims, uvw);
    }


    public static double sigmoidFun(double rad, double distance, double gamma) {
        double x = -gamma * (distance - rad);
        return Math.exp(x) / (1. + Math.exp(x));
    }

    /**
     * Set all vertical speeds in the grid to ensure zero net balance in all
     * cells. Assume no flow through the top of the uppermost layer, and free
     * flow through the bottom of the lowermost layer.
     * @param dims Domain dimensions
     * @param uvw The full current field
     */
    public static void computeVerticalSpeeds(int[] dims, double[][][][] uvw) {

        for (int i = 0; i < dims[0]; i++)
            for (int j = 0; j < dims[1]; j++)
                for (int k = 0; k < dims[2]; k++) {
                    // Balance between horizontal flows in and out of this cell:
                    double balance = uvw[i][j][k][0] - uvw[i + 1][j][k][0]
                            + uvw[i][j][k][1] - uvw[i][j + 1][k][1];
                    double wIn = 0;
                    // If we are not in the top layer we may have an inflow from above:
                    wIn += uvw[i][j][k][2];

                    // Set downwards outflow to compensate for imbalance:
                    uvw[i][j][k+1][2] = balance + wIn;
                }

    }


    public static void stats(double[][][][] hydro) {
        double maxvalX = Double.MIN_VALUE, minvalX = Double.MAX_VALUE;
        double maxvalY = Double.MIN_VALUE, minvalY = Double.MAX_VALUE;
        double maxvalZ = Double.MIN_VALUE, minvalZ = Double.MAX_VALUE;
        double maxBal = Double.MIN_VALUE, minBal = Double.MAX_VALUE;
        int cellCount = 0, imbalancedCount = 0;

        for (int i=0; i<hydro.length; i++)
            for (int j=0; j<hydro[i].length; j++)
                for (int k=0; k<hydro[i][j].length; k++) {

                    // Balance:
                    if ((i < hydro.length-1) && (j < hydro[i].length-1) && (k < hydro[i][j].length-1)) {
                        double balance = hydro[i][j][k][0] - hydro[i + 1][j][k][0]
                                + hydro[i][j][k][1] - hydro[i][j + 1][k][1]
                                + hydro[i][j][k][2] - hydro[i][j][k + 1][2];
                        if (balance > maxBal)
                            maxBal = balance;
                        if (balance < minBal)
                            minBal = balance;

                        cellCount++;
                        if (Math.abs(balance) > 1e-11)
                            imbalancedCount++;
                    }

                    // Min and max values:
                    if (hydro[i][j][k][0] > maxvalX)
                        maxvalX = hydro[i][j][k][0];
                    if (hydro[i][j][k][1] > maxvalY)
                        maxvalY = hydro[i][j][k][1];
                    if (hydro[i][j][k][2] > maxvalZ)
                        maxvalZ = hydro[i][j][k][2];
                    if (hydro[i][j][k][0] < minvalX)
                        minvalX = hydro[i][j][k][0];
                    if (hydro[i][j][k][1] < minvalY)
                        minvalY = hydro[i][j][k][1];
                    if (hydro[i][j][k][2] < minvalZ)
                        minvalZ = hydro[i][j][k][2];
                }
        System.out.println("U values: "+minvalX +" - "+maxvalX);
        System.out.println("V values: "+minvalY +" - "+maxvalY);
        System.out.println("W values: "+minvalZ +" - "+maxvalZ);
        System.out.println("Balance: "+minBal +" - "+maxBal);
        System.out.println("Imbalanced / total: "+imbalancedCount + " / "+cellCount);
    }

    public static void main(String[] args) {
        double dx = 2;
        double rad = 6;
        double currRedFactor = 1. - 0.8;
        int[] dims = new int[] {12, 12, 5};
        double[][] profile = new double[dims[2]][2];
        double[][][][] uvw = new double[dims[0]+1][dims[1]+1][dims[2]+1][3];

        for (int k=0; k<dims[2]; k++) {
            profile[k][0] = 0.08;
            profile[k][1] = 0.;
        }
        getBjoroyaHydraulicField(dims, dx, rad, currRedFactor, profile, uvw);

        NumberFormat nf = NumberFormat.getNumberInstance(Locale.ENGLISH);
        nf.setMaximumFractionDigits(4);
        nf.setMinimumFractionDigits(4);

        System.out.println("U field (k=0):");
        for (int i=0; i<dims[0]+1; i++) {
            for (int j = 0; j < dims[1]; j++) {
                System.out.print(nf.format(uvw[i][j][0][0]) + "  ");
            }
            System.out.println();
        }

        System.out.println("V field (k=0):");
        for (int i=0; i<dims[0]; i++) {
            for (int j = 0; j < dims[1]+1; j++) {
                System.out.print(nf.format(uvw[i][j][0][1]) + "  ");
            }
            System.out.println();
        }

        System.out.println("W field (k=4):");
        for (int i=0; i<dims[0]; i++) {
            for (int j = 0; j < dims[1]; j++) {
                System.out.print(nf.format(uvw[i][j][3][2]) + "  ");
            }
            System.out.println();
        }
    }
}
