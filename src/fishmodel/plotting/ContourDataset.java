package fishmodel.plotting;

import org.jfree.data.xy.AbstractXYZDataset;

/**
 * Created by malv on 04.05.15.
 */
public class ContourDataset extends AbstractXYZDataset {

    private double[][][] pConc = null;
    private int sliceDim = 0;
    private int sliceCoord = 0;
    private int[] pdim = null, adim = null;
    private boolean yReverse = false;

    public ContourDataset() {

    }

    public void setData(double[][][] pConc, int sliceDim, int sliceCoord) {
        this.pConc = pConc;
        this.sliceCoord = sliceCoord;
        this.sliceDim = sliceDim;
        yReverse = (sliceDim < 2);
        int[] dims = new int[] {pConc.length, pConc[0].length, pConc[0][0].length};
        int piv=0;
        adim = new int[2];
        pdim = new int[2];
        for (int i = 0; i < dims.length; i++) {
            if (i != sliceDim) {
                adim[piv] = i;
                pdim[piv] = dims[i];
                piv++;
            }
        }

        /*for (int i = 0; i < adim.length; i++) {

            System.out.println("Adim["+i+"] = "+adim[i]);
            System.out.println("Pdim["+i+"] = "+pdim[i]);
        } */
    }

    @Override
    public int getSeriesCount() {
        return 1;
    }

    @Override
    public Comparable getSeriesKey(int series) {
        return String.valueOf(series);
    }

    public int getMaxX() {
        if (pdim != null)
            return pdim[0];
        else return 0;
    }

    public int getMaxY() {
        if (pdim != null)
            return pdim[1];
        else return 0;
    }

    @Override
    public int getItemCount(int series) {
        if (pdim != null)
            return pdim[0]*pdim[1];
        else return 0;
    }

    @Override
    public Number getZ(int series, int item) {

        int x = getX(series, item).intValue();
        int y = getYInternal(item, false).intValue();
        if (pdim != null) {
            if (sliceDim == 0) {
                return pConc[sliceCoord][x][y];
            }
            else if (sliceDim == 1) {
                return pConc[x][sliceCoord][y];
            }
            else {
                return pConc[x][y][sliceCoord];
            }
        }
        else
            return null;
    }


    @Override
    public Number getX(int series, int item) {
        if (pdim != null) {
            return item % pdim[0];
        }
        else
            return null;
    }

    @Override
    public Number getY(int series, int item) {
        return getYInternal(item, yReverse);
    }

    public Number getYInternal(int item, boolean yReverse) {
        if (pdim != null) {
            return yReverse ? pdim[1]-1-item/pdim[0] : item/pdim[0];
        }
        else
            return null;
    }
}
