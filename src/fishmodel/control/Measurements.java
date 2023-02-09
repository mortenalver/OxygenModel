package fishmodel.control;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Random;

/**
 * Created by malv on 09.06.17.
 */
public class Measurements {

    private boolean noiseActivated = true;
    private int[] cageDims;
    private double dxy;
    private double[] feedByDepth;
    private double[][][] feedConcentration;
    private Random rnd;
    private File outFile = new File("C:\\Temp\\meas.csv");
    static boolean first = true;

    public Measurements(int[] cageDims, double dxy, double[] feedByDepth, double[][][] feedConcentration) {
        this.cageDims = cageDims;
        this.dxy = dxy;
        this.feedByDepth = feedByDepth;
        this.feedConcentration = feedConcentration;
        rnd = new Random();
    }

    private void logDataPoint(double precise, double noisy) {
        try {
            FileWriter out = null;
            if (first) {
                out = new FileWriter(outFile);
                first = false;
            }
            else
                out = new FileWriter(outFile, true);

            out.write(precise+"\t"+noisy+"\r\n");

            try { out.close(); } catch (IOException ex) { ex.printStackTrace(); }
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public void setNoiseActivated(boolean noiseActivated) {
        this.noiseActivated = noiseActivated;
    }

    public double measurePelletConcentrationAtPoint(int[] samplePoint) {
        double value = feedConcentration[samplePoint[0]][samplePoint[1]][samplePoint[2]] /(dxy*dxy*dxy);
        double precise = value;
        if (noiseActivated) {
            double noise = rnd.nextGaussian()*Math.max(Math.min(value/2.5, 0.4), 0.025);
            value += noise;
            value = Math.max(0, value);
        }
        //logDataPoint(precise, value);
        return value;
    }

    public double measureTotalFeed() {
        double value = measureTotalFeedExact();
        double precise = value;
        if (noiseActivated) {
            double std = Math.max(Math.min(value/2.5, 3000), 100);//Math.min(value, 200);
            double noise = rnd.nextGaussian()*std;
            value += noise;
            value = Math.max(0, value);
        }
        //logDataPoint(precise, value);
        return value;
    }

    private double measureTotalFeedExact() {
        double totFeed = 0;
        for (int i=0; i<cageDims[2]; i++) {
            totFeed += feedByDepth[i];
        }
        return totFeed;
    }

    public double measureAverageFeedDepth() {
        double totFeed = measureTotalFeedExact();
        double value = 0;
        for (int i=0; i<cageDims[2]; i++) {
            value += feedByDepth[i]*(((double)i)+0.5)*dxy;
        }
        if (totFeed > 0)
            value /= totFeed;
        else value = dxy;
        double precise = value;
        if (noiseActivated) {
            double std = 1.5;
            double noise = rnd.nextGaussian()*std;
            value += noise;
            value = Math.max(0, value);
        }
        logDataPoint(precise, value);
        return value;
    }
}
