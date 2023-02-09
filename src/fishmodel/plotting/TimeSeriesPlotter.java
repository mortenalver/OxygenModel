package fishmodel.plotting;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.data.time.TimeSeriesCollection;

import java.awt.*;

/**
 * Created by malv on 07.05.15.
 */
public class TimeSeriesPlotter {

    ChartPanel panel;

    public TimeSeriesPlotter(TimeSeriesCollection ts, String title, String xLabel, String yLabel, boolean inverted) {
        this.panel = initChart(ts, title, xLabel, yLabel, inverted);

    }

    public ChartPanel getPanel() {
        return panel;
    }


    private ChartPanel initChart(TimeSeriesCollection ts, String title, String xLabel, String yLabel, boolean inverted) {

        JFreeChart fc = ChartFactory.createTimeSeriesChart
                (title, xLabel, yLabel, ts, true, true, false);
        //fc.removeLegend();
        ChartPanel panel = new ChartPanel(fc);
        if (inverted)
            fc.getXYPlot().getRangeAxis().setInverted(true);
        Font tickLabelFont = new Font("plain", Font.PLAIN, 12);
        fc.getXYPlot().getRangeAxis().setTickLabelFont(tickLabelFont);
        //fc.getLegend().setVisible(false);
        // Increase the maximum draw width and height for the chart panel in order
        // to prevent scaling of the fonts:
        panel.setMaximumDrawWidth(2000);
        panel.setMaximumDrawHeight(2000);

        fc.getXYPlot().getDomainAxis().setTickLabelFont(tickLabelFont);
        ValueAxis xAxis = fc.getXYPlot().getDomainAxis();
        xAxis.setAutoRange(true);

        ValueAxis yAxis = fc.getXYPlot().getRangeAxis();
        //yAxis.setAutoRangeMinimumSize(1.0);
        panel.setPreferredSize(new Dimension(1100, 500));
        return panel;
    }

}
