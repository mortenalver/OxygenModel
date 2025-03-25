package fishmodel;

import java.util.ArrayList;
import java.util.Arrays;

public class MultiCageUtils {

    protected static double[] binEdges = null;

    public static double[] getBinEdges() {
        if (binEdges == null) {
            // Define bin edges first time:
            double min = 0., step = 0.1, max = 10.;
            int steps = (int)Math.round((max-min)/step);
            binEdges = new double[steps];
            for (int i = 0; i < binEdges.length; i++) {
                binEdges[i] = min + i*step;
            }
        }
        return binEdges;
    }

    public static int[][] getRanges(double[] pos, double dxy, double rad) {
        int[][] ranges = new int[2][2];
        int radI = (int)(Math.ceil(rad/dxy));
        ranges[0] = new int[] {(int)Math.floor(pos[0]/dxy)-radI-1,
                (int)Math.floor(pos[0]/dxy)+radI+1};
        ranges[1] = new int[] {(int)Math.floor(pos[1]/dxy)-radI-1,
                (int)Math.floor(pos[1]/dxy)+radI+1};
        return ranges;
    }

    public static double[] getO2AffSums(ArrayList<double[]> cagePos, double dxy, double rad, double[][][] o2Affinity) {
        double[] affSums = new double[cagePos.size()];
        for (int cg=0; cg<affSums.length; cg++) {
            int[][] ranges = getRanges(cagePos.get(cg), dxy, rad);
            int x0 = ranges[0][0],
                x1 = ranges[0][1],
                y0 = ranges[1][0],
                y1 = ranges[1][1];
            for (int i=x0; i<x1; i++)
                for (int j=y0; j<y1; j++)
                    for (int k=0; k<o2Affinity[0][0].length; k++)
                        affSums[cg] += o2Affinity[i][j][k];
        }
        return affSums;
    }

    public static ArrayList<CageStats> getCageStats(double[][][] o2, boolean[][][] mask, ArrayList<double[]> cagePos, double cageRad,
                                          double dxy, double hypoxiaThreshold) {
        int radI = (int)(Math.ceil(cageRad/dxy));

        ArrayList<CageStats> res = new ArrayList<CageStats>(cagePos.size());

        for (int i = 0; i < cagePos.size(); i++) {
            double[] pos = cagePos.get(i);
            int[][] ranges = getRanges(pos, dxy, cageRad);
            int[] xrange = ranges[0];
            int[] yrange = ranges[1];
            /*int[] xrange = new int[] {(int)Math.floor(pos[0]/dxy)-radI-1,
                    (int)Math.floor(pos[0]/dxy)+radI+1};
            int[] yrange = new int[] {(int)Math.floor(pos[1]/dxy)-radI-1,
                    (int)Math.floor(pos[1]/dxy)+radI+1};*/
            CageStats cStats = getStats(o2, mask, xrange, yrange, hypoxiaThreshold);
            res.add(cStats);


        }
        return res;
    }

    public static ArrayList<CageStats> getCageStatsRectangular(double[][][] o2, boolean[][][] mask, ArrayList<double[]> cagePos, double[] unitSizeM,
                                                    double dxy, double hypoxiaThreshold) {
        int[] unitSizeI = new int[] {(int)(Math.ceil(unitSizeM[0]/dxy)), (int)(Math.ceil(unitSizeM[1]/dxy))};
        //System.out.println("unit size gridpoints: "+unitSizeI[0]+" / "+unitSizeI[1]);
        ArrayList<CageStats> res = new ArrayList<CageStats>(cagePos.size());

        for (int i = 0; i < cagePos.size(); i++) {
            double[] pos = cagePos.get(i);
            //System.out.println("Cage pos : "+pos[0]+" / "+pos[1]);
            int[] xrange = new int[] {(int)Math.floor(pos[0]/dxy)-unitSizeI[0]/2-1,
                    (int)Math.floor(pos[0]/dxy)+unitSizeI[0]/2+1};
            int[] yrange = new int[] {(int)Math.floor(pos[1]/dxy)-unitSizeI[1]/2-1,
                    (int)Math.floor(pos[1]/dxy)+unitSizeI[1]/2+1};
            //System.out.println("xrange: "+xrange[0]+" - "+xrange[1]);
            //System.out.println("yrange: "+yrange[0]+" - "+yrange[1]);
            CageStats cStats = getStats(o2, mask, xrange, yrange, hypoxiaThreshold);
            res.add(cStats);

        }
        return res;
    }

    // Get min, 5% percentile, 10% percentile, max, mean and hypoxic fraction for cage
    protected static CageStats getStats(double[][][] o2, boolean[][][] mask, int[] xrange, int[] yrange,
                                              double hypoxiaThreshold) {

        int totCells = 0;
        int cellsBelow = 0;
        double max = Double.MIN_VALUE,
                min = Double.MAX_VALUE,
                sum = 0;
        double[] allVals = new double[(1+xrange[1]-xrange[0])*(1+yrange[1]-yrange[0])*o2[0][0].length];
        //System.out.println("xrange: "+xrange[0]+" - "+xrange[1]);
        //System.out.println("yrange: "+yrange[0]+" - "+yrange[1]);
        //System.out.println("Allvals size: "+allVals.length);
        for (int i=xrange[0]; i<=xrange[1]; i++)
            for (int j=yrange[0]; j<=yrange[1]; j++)
                for (int k=0; k<o2[i][j].length; k++) {
                    if (mask[i][j][k]) {
                        allVals[totCells] = o2[i][j][k];
                        totCells++;
                        if (o2[i][j][k] < hypoxiaThreshold)
                            cellsBelow++;
                        if (o2[i][j][k] > max)
                            max = o2[i][j][k];
                        if (o2[i][j][k] < min)
                            min = o2[i][j][k];
                        sum += o2[i][j][k];
                    }
                }

        double[] allVals2 = new double[totCells];
        System.arraycopy(allVals, 0, allVals2, 0, totCells);
        Arrays.sort(allVals2);

        //long tic = System.currentTimeMillis();
        int[] histValues = getHistValues(allVals2);

        int sumBins = 0;
        for (int i = 0; i < histValues.length; i++) {
            sumBins += histValues[i];
            //System.out.println(">> bin "+i+": "+histValues[i]);
        }
        //long toc = System.currentTimeMillis();

        //System.out.println("Sum bin values: "+sumBins);
        //System.out.println("Number of values: "+allVals2.length);

        return new CageStats(new double[] {min, allVals2[totCells/20], allVals2[totCells/10], max, sum/((double)totCells), ((double)cellsBelow)/((double)totCells)},
                histValues);
    }

    protected static int[] getHistValues(double[] allVals) {

        double[] binEdg = getBinEdges();
        int[] histVals = new int[binEdg.length];

        // Iterate over all values and find how many go into each bin progressively:
        int piv = 0;
        int counter = 0;
        for (int i = 0; i < binEdg.length-1; i++) {
            // Count until we reach an element that doesn't belong in this bin:
            while ((piv < allVals.length) && (allVals[piv] < binEdg[i+1])) {
                histVals[i]++;
                piv++;
            }

        }
        // Put the rest of the elements in the last bin:
        histVals[binEdg.length-1] = allVals.length-piv;



        return histVals;
    }
}
