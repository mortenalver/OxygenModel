package fishmodel.plotting;

import fishmodel.CageMasking;
import fishmodel.pellets.PelletSpreaderModel;
import org.jfree.chart.ChartMouseEvent;
import org.jfree.chart.ChartMouseListener;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.XYShapeAnnotation;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.GrayPaintScale;
import org.jfree.chart.renderer.xy.XYBlockRenderer;

import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;

/**
 * Created by malv on 21.05.15.
 */
public class FeedSettingsPanel {

    protected JPanel panel;
    protected NumberAxis valueAxis, domainAxis;
    XYBlockRenderer blocRenderer;
    XYPlot xyp;
    JFreeChart chart;
    ChartPanel chartPanel;
    JSlider airSpeedS, feedingRateS, currentSpdS, currentDirS;

    public static final int
        AIR_SPEED_MIN = 10,
        AIR_SPEED_MAX = 30,
        FEEDING_RATE_MAX = 100,
        CURRENTSPEED_MAX = 100;

    protected int centerX, centerY;

    protected double[][] feedingRate;
    protected boolean[][][] mask;
    protected SurfaceDataset dataset = null;
    private int[] cageDims;
    private double rad;
    private double dxy;
    private int spreaderType;
    private boolean tiltUp;
    private int airSpeed;
    private double feedingRateMult;
    private double[] windSpeed;

    public FeedSettingsPanel(int[] cageDims, int centerX, int centerY, double rad, double dxy, int spreaderType,
                             boolean tiltUp, int airSpeed, double feedingRateMult, double[] windSpeed, double[] currentSpeed) {
        this.cageDims = cageDims;
        this.rad = rad;
        this.dxy = dxy;
        this.spreaderType = spreaderType;
        this.tiltUp = tiltUp;
        this.airSpeed = airSpeed;
        this.feedingRateMult = feedingRateMult;
        this.windSpeed = new double[windSpeed.length];
        for (int i = 0; i < windSpeed.length; i++) {
            this.windSpeed[i] = windSpeed[i];

        }


        panel = new JPanel();
        panel.setLayout(new BorderLayout());
        this.centerX = centerX;
        this.centerY = centerY;

        blocRenderer = new XYBlockRenderer();
        blocRenderer.setPaintScale(new GrayPaintScale(0, 0.5));

        mask = CageMasking.circularMasking(cageDims, dxy, rad, false);
        feedingRate = new double[cageDims[0]][cageDims[1]];

        calculateFeedingRate();

        dataset = new SurfaceDataset(feedingRate, false);
        valueAxis = new NumberAxis("En en en");
        domainAxis = new NumberAxis("To to to");
        domainAxis.setVisible(false);
        valueAxis.setVisible(false);

        xyp = new XYPlot(dataset, domainAxis, valueAxis, blocRenderer);

        // Add annotation showing the cage boundary:
        XYShapeAnnotation a1 = new XYShapeAnnotation(
                new Ellipse2D.Double((cageDims[0]/2)-(rad/dxy)-0.5, (cageDims[1]/2)-(rad/dxy)-0.5, 2*rad/dxy, 2*rad/dxy), new BasicStroke(2.0f), Color.white);
        xyp.addAnnotation(a1);
        chart = new JFreeChart(xyp);
        chart.getLegend().setVisible(false);
        chartPanel = new ChartPanel(chart);
        chartPanel.addChartMouseListener(new ClickListener());
        chartPanel.setPreferredSize(new Dimension(500,500));
        JLabel asLab = new JLabel("Air speed (m/s)");
        airSpeedS = new JSlider(AIR_SPEED_MIN, AIR_SPEED_MAX, airSpeed);
        airSpeedS.setOrientation(JSlider.VERTICAL);
        airSpeedS.addChangeListener(new AirSpeedChangeListener());
        //Turn on labels at major tick marks.
        airSpeedS.setMajorTickSpacing(5);
        airSpeedS.setMinorTickSpacing(1);
        airSpeedS.setPaintTicks(true);
        airSpeedS.setPaintLabels(true);
        JPanel epan = new JPanel();
        epan.setLayout(new BorderLayout());
        epan.add(asLab, BorderLayout.NORTH);
        epan.add(airSpeedS, BorderLayout.CENTER);

        JLabel frLab = new JLabel("Feeding rate (kg/min)");
        feedingRateS = new JSlider(0, FEEDING_RATE_MAX, (int)feedingRateMult);
        feedingRateS.setOrientation(JSlider.VERTICAL);
        feedingRateS.addChangeListener(new FeedingRateChangeListener());
        //Turn on labels at major tick marks.
        feedingRateS.setMajorTickSpacing(50);
        feedingRateS.setMinorTickSpacing(10);
        feedingRateS.setPaintTicks(true);
        feedingRateS.setPaintLabels(true);
        JPanel fpan = new JPanel();
        fpan.setLayout(new BorderLayout());
        fpan.add(frLab, BorderLayout.NORTH);
        fpan.add(feedingRateS, BorderLayout.CENTER);

        JLabel csLab = new JLabel("Current speed (cm/s)");
        currentSpdS = new JSlider(0, CURRENTSPEED_MAX,
                (int)(100*Math.sqrt(currentSpeed[0]*currentSpeed[0] + currentSpeed[1]*currentSpeed[1])));
        currentSpdS.setOrientation(JSlider.VERTICAL);
        //Turn on labels at major tick marks.
        currentSpdS.setMajorTickSpacing(20);
        currentSpdS.setMinorTickSpacing(5);
        currentSpdS.setPaintTicks(true);
        currentSpdS.setPaintLabels(true);
        JPanel cspan = new JPanel();
        cspan.setLayout(new BorderLayout());
        cspan.add(csLab, BorderLayout.NORTH);
        cspan.add(currentSpdS, BorderLayout.CENTER);

        JLabel cdLab = new JLabel("Current direction");
        int oldCD = (int)((180./Math.PI)*Math.atan2(currentSpeed[1], currentSpeed[0]));
        oldCD = 90-oldCD;
        if (oldCD < 0) oldCD += 360;

        System.out.println("OldCD: "+oldCD);
        currentDirS = new JSlider(0, 360, oldCD);
        currentDirS.setOrientation(JSlider.VERTICAL);
        //Turn on labels at major tick marks.
        currentDirS.setMajorTickSpacing(90);
        currentDirS.setMinorTickSpacing(30);
        currentDirS.setPaintTicks(true);
        currentDirS.setPaintLabels(true);
        JPanel cdpan = new JPanel();
        cdpan.setLayout(new BorderLayout());
        cdpan.add(cdLab, BorderLayout.NORTH);
        cdpan.add(currentDirS, BorderLayout.CENTER);

        JPanel currentPan = new JPanel();
        currentPan.setLayout(new BorderLayout());
        currentPan.add(cspan, BorderLayout.WEST);
        currentPan.add(cdpan, BorderLayout.CENTER);

        JPanel slipan = new JPanel();
        slipan.setLayout(new BorderLayout());
        slipan.add(epan, BorderLayout.WEST);
        slipan.add(fpan, BorderLayout.CENTER);
        slipan.add(currentPan, BorderLayout.EAST);

        panel.add(slipan, BorderLayout.EAST);
        panel.add(chartPanel, BorderLayout.CENTER);

        /*Plot currentDirPlot = CurrentDirPlot.getPolarPlot(this.currentSpeed);
        JFreeChart cdChart = new JFreeChart(currentDirPlot);
        ChartPanel ccp = new ChartPanel(cdChart);
        ccp.addChartMouseListener(new CurrentSpeedClickListener());
        panel.add(ccp, BorderLayout.WEST);*/

    }

    public void calculateFeedingRate() {
        PelletSpreaderModel.setPelletDist(feedingRate, centerX, centerY, dxy, 0, airSpeed, spreaderType, tiltUp, windSpeed, mask);
        setColorScaleBounds(0, maxVal(feedingRate));
        if (dataset != null)
            dataset.setData(feedingRate);
    }



    public void setColorScaleBounds(double lower, double upper) {
        blocRenderer.setPaintScale(new RedBlueScale(lower, upper));
    }

    private double maxVal(double[][] matr) {
        double max = 0;
        for (int i = 0; i < matr.length; i++) {
            for (int j = 0; j < matr[i].length; j++) {
                if (matr[i][j] > max)
                    max = matr[i][j];
            }
        }
        return max;
    }
    
    public JPanel getPanel() {
        return panel;
    }

    class ClickListener implements ChartMouseListener {

        @Override
        public void chartMouseClicked(ChartMouseEvent e) {
            Point2D p = chartPanel.translateScreenToJava2D(e.getTrigger().getPoint());
            Rectangle2D plotArea = chartPanel.getScreenDataArea();
            XYPlot plot = (XYPlot) chart.getPlot(); // your plot
            double chartX = plot.getDomainAxis().java2DToValue(p.getX(), plotArea, plot.getDomainAxisEdge());
            double chartY = plot.getRangeAxis().java2DToValue(p.getY(), plotArea, plot.getRangeAxisEdge());
            //System.out.println("chartX = "+chartX+", chartY = "+chartY);
            final int cx = (int)(Math.round(chartX)), cy = (int)(Math.round(chartY));
            // Update feed distribution on a separate thread:
            new Thread(new Runnable() {
                @Override
                public void run() {
                    centerX = cx;
                    centerY = cy;
                    calculateFeedingRate();
                }
            }).start();
        }

        @Override
        public void chartMouseMoved(ChartMouseEvent e) {

        }

    }



    class AirSpeedChangeListener implements ChangeListener {
        @Override
        public void stateChanged(ChangeEvent e) {
            airSpeed = airSpeedS.getValue();
            new Thread(new Runnable() {
                @Override
                public void run() {
                    calculateFeedingRate();
                }
            }).start();

        }
    }

    class FeedingRateChangeListener implements ChangeListener {
        @Override
        public void stateChanged(ChangeEvent e) {
            feedingRateMult = feedingRateS.getValue();
        }
    }
    protected static boolean okPressed = false;

    public boolean showFeedSettingsDialog() {

        okPressed = false;

        final JDialog diag = new JDialog((JFrame)null, "Feeding setup", true);
        diag.getContentPane().setLayout(new BorderLayout());
        diag.getContentPane().add(getPanel(), BorderLayout.CENTER);

        JPanel butPan = new JPanel();
        JButton ok = new JButton("Ok"), cancel = new JButton("Cancel");
        butPan.add(ok); butPan.add(cancel);

        ok.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                okPressed = true;
                diag.dispose();
            }
        });
        cancel.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                okPressed = false;
                diag.dispose();
            }
        });

        diag.getContentPane().add(butPan, BorderLayout.SOUTH);
        diag.pack();
        diag.setLocationRelativeTo(null);
        diag.setVisible(true);

        return okPressed;
    }

    public int getCenterX() {
        return centerX;
    }

    public int getCenterY() {
        return centerY;
    }

    public int getAirSpeed() {
        return airSpeed;
    }

    public double getFeedingRateMult() { return feedingRateMult; }

    public double[] getCurrentSpeed() {
        double[] result = new double[2];
        double speedVal = 0.01*(double)(currentSpdS.getValue());
        result[0] = speedVal*Math.sin((Math.PI/180.)*(double)currentDirS.getValue());
        result[1] = speedVal*Math.cos((Math.PI/180.)*(double)currentDirS.getValue());
        return result;

    }

}
