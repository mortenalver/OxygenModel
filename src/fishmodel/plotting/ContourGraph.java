package fishmodel.plotting;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.XYShapeAnnotation;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.GrayPaintScale;
import org.jfree.chart.renderer.xy.XYBlockRenderer;

import java.awt.*;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Line2D;

/**
 * Created by malv on 04.05.15.
 */
public class ContourGraph {

    ContourDataset dataset;
    XYPlot xyp;
    XYBlockRenderer blocRenderer;
    JFreeChart chart;
    ChartPanel panel;
    ValueAxis valueAxis, domainAxis;
    private int sliceDim;
    private int sliceCoord;

    private XYShapeAnnotation hotSpotAnnotation = null;

    public ContourGraph() {
        dataset = new ContourDataset();
        valueAxis = new NumberAxis("En en en");
        domainAxis = new NumberAxis("To to to");
        domainAxis.setVisible(false);
        valueAxis.setVisible(false);
        blocRenderer = new XYBlockRenderer();
        blocRenderer.setPaintScale(new GrayPaintScale(0, 50));
        xyp = new XYPlot(dataset, domainAxis, valueAxis, blocRenderer);

        chart = new JFreeChart(xyp);
        panel = new ChartPanel(chart);
        chart.getLegend().setVisible(false);
    }

    public void addCageAnnotation(int centerX, int centerY, double rad, String label) {
        // Add annotation showing the cage boundary from above:
        XYShapeAnnotation a1 = new XYShapeAnnotation(
                new Ellipse2D.Double(centerX-rad-0.5, centerY-rad-0.5, 2*rad, 2*rad), new BasicStroke(2.0f), Color.white);
        xyp.addAnnotation(a1);

        XYTextAnnotation t1 = new XYTextAnnotation(label, centerX-rad+0.85, centerX-rad-0.75);
        t1.setFont(new Font("plain", Font.PLAIN, 24));
        t1.setPaint(Color.white);
        xyp.addAnnotation(t1);
    }

    public void addSidewaysCageAnnotation(int centerX, double rad, double cylDepth, double totDepth, double modelBottom) {
        // Add annotation showing the cage silhouette from the side:
        System.out.println(centerX-rad-0.5);
        System.out.println("CylDepth = "+cylDepth);
        XYShapeAnnotation a1 = new XYShapeAnnotation(
                new Line2D.Double(centerX-rad-0.5, modelBottom, centerX-rad-0.5+0.0001, modelBottom-cylDepth), new BasicStroke(2.0f), Color.white);
        xyp.addAnnotation(a1);
        a1 = new XYShapeAnnotation(
                new Line2D.Double(centerX+rad-0.5, modelBottom, centerX+rad-0.5+0.0001, modelBottom-cylDepth), new BasicStroke(2.0f), Color.white);
        xyp.addAnnotation(a1);
        a1 = new XYShapeAnnotation(
                new Line2D.Double(centerX+rad-0.5, modelBottom-cylDepth, centerX-0.5, modelBottom-totDepth), new BasicStroke(2.0f), Color.white);
        xyp.addAnnotation(a1);
        a1 = new XYShapeAnnotation(
                new Line2D.Double(centerX-rad-0.5, modelBottom-cylDepth, centerX-0.5, modelBottom-totDepth), new BasicStroke(2.0f), Color.white);
        xyp.addAnnotation(a1);

        //new Ellipse2D.Double(centerX-rad-0.5, centerY-rad-0.5, 2*rad, 2*rad), new BasicStroke(2.0f), Color.white);
    }

    public void setHotSpotAnnotation(double x, double y, double z, double modelBottom) {
        if (hotSpotAnnotation != null)
            xyp.removeAnnotation(hotSpotAnnotation);
        double size = 0.2;
        if (sliceDim == 0) {
            hotSpotAnnotation = new XYShapeAnnotation(new Ellipse2D.Double(y-size,modelBottom-z-size, 2*size, 2*size));
        } else if (sliceDim == 1) {
            hotSpotAnnotation = new XYShapeAnnotation(new Ellipse2D.Double(x-size,modelBottom-z-size, 2*size, 2*size));
        } else {
            hotSpotAnnotation = new XYShapeAnnotation(new Ellipse2D.Double(x-size,y-size, 2*size, 2*size));
        }

        xyp.addAnnotation(hotSpotAnnotation);
    }

    public void setColorScaleBounds(double lower, double upper) {
        //blocRenderer.setPaintScale(new GrayPaintScale(lower, upper));
        blocRenderer.setPaintScale(new RedBlueScale(lower, upper));
    }

    public void dataUpdated() {
        chart.fireChartChanged();
    }

    public void setSlice(int sliceDim, int sliceCoord) {

        this.sliceDim = sliceDim;
        this.sliceCoord = sliceCoord;
    }

    public void setData(final double[][][] pConc) {
        dataset.setData(pConc, sliceDim, sliceCoord);
        domainAxis.setRange(0, dataset.getMaxX()-0.5);
        valueAxis.setRange(0, dataset.getMaxY()-0.5);
        chart.fireChartChanged();
    }

    public ChartPanel getPanel() {
        return panel;
    }
}
