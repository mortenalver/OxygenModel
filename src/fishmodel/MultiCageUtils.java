package fishmodel;

import java.util.ArrayList;
import java.util.Arrays;

public class MultiCageUtils {

    // Get min, 5% percentile, 10% percentile, max, mean and hypoxic fraction for cage
    protected static double[] getStats(double[][][] o2, boolean[][][] mask, int[] xrange, int[] yrange,
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

        return new double[] {min, allVals2[totCells/20], allVals2[totCells/10], max, sum/((double)totCells), ((double)cellsBelow)/((double)totCells)};
    }
    public static double[][] getCageStats(double[][][] o2, boolean[][][] mask, ArrayList<double[]> cagePos, double cageRad,
                                        double dxy, double hypoxiaThreshold) {
        double[][] res = new double[cagePos.size()][4];
        int radI = (int)(Math.ceil(cageRad/dxy));

        for (int i = 0; i < res.length; i++) {
            double[] pos = cagePos.get(i);
            int[] xrange = new int[] {(int)Math.floor(pos[0]/dxy)-radI-1,
                    (int)Math.floor(pos[0]/dxy)+radI+1};
            int[] yrange = new int[] {(int)Math.floor(pos[1]/dxy)-radI-1,
                    (int)Math.floor(pos[1]/dxy)+radI+1};
            double[] stats = getStats(o2, mask, xrange, yrange, hypoxiaThreshold);
            res[i] = stats;


        }
        return res;
    }
}
