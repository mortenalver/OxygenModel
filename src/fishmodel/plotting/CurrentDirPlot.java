package fishmodel.plotting;

import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.PolarPlot;
import org.jfree.chart.renderer.DefaultPolarItemRenderer;
import org.jfree.data.xy.DefaultXYDataset;

/**
 * Created by malv on 09.06.15.
 */
public class CurrentDirPlot {

    public static final double MAX_SPEED = 2.0;

    public static PolarPlot getPolarPlot(double[] currentVector) {

        DefaultXYDataset dataset = new DefaultXYDataset();
        NumberAxis valueAxis = new NumberAxis("Current speed");
        valueAxis.setAutoRange(false);
        valueAxis.setRange(0, MAX_SPEED);
        double angle = Math.atan2(currentVector[1], currentVector[0]);
        double speed = Math.sqrt(currentVector[0]*currentVector[0] + currentVector[1]*currentVector[1]);

        double[][] data = new double[2][1];
        data[0][0] = angle*180/Math.PI;
        data[1][0] = speed;
        dataset.addSeries("Test", data);
        return new PolarPlot(dataset, valueAxis, new DefaultPolarItemRenderer());


    }

}
