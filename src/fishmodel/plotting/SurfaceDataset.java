package fishmodel.plotting;

import org.jfree.data.xy.AbstractXYZDataset;

/**
 * Created by malv on 04.05.15.
 */
public class SurfaceDataset extends AbstractXYZDataset {

    private double[][] pConc;
    private int maxX=0, maxY=0, nItems=0;
    private boolean yReverse;

    public SurfaceDataset(double[][] pConc, boolean yReverse) {
        this.yReverse = yReverse;
        setData(pConc);
    }

    public void setData(double[][] pConc) {
        this.pConc = pConc;
        maxX = pConc.length;
        maxY = pConc[0].length;
        nItems = maxX*maxY;
        fireDatasetChanged();
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
        return maxX;
    }

    public int getMaxY() {
        return maxY;
    }

    @Override
    public int getItemCount(int series) {
        return nItems;
    }

    @Override
    public Number getZ(int series, int item) {

        int x = getX(series, item).intValue();
        int y = getYInternal(item, false).intValue();
        return pConc[x][y];
    }


    @Override
    public Number getX(int series, int item) {
        return item % maxX;

    }

    @Override
    public Number getY(int series, int item) {
        return getYInternal(item, yReverse);
    }

    public Number getYInternal(int item, boolean yReverse) {
        return yReverse ? maxY-1-item/maxX : item/maxX;

    }
}
