package fishmodel.control;

/**
*
 */
public class FeedRateController {

    // Initial feeding rate (kg/min):
    public static double nominalFeedRate = 28;
    /**
     * This is the method that computes the feeding rate.
     * @param currentFeedRate The current feeding rate
     * @param cageDims Grid dimensions (horizontal 1, horizontal 2 and vertical)
     * @param dxy Grid resolution (m)
     * @param meas Object providing access to measurements
     * @return The updated feeding rate
     */
    public double getFeedRate(double currentFeedRate, int[] cageDims, double dxy, Measurements meas) {

        // Example calculation of the depth of the lower edge of the grid:
        //double maxDepth = ((double)cageDims[2])*dxy;

        double measurement = 0;

        // -------------------------------------------------------------------------------------------------------------
        // Define the measurement available to your controller here:
        // -------------------------------------------------------------------------------------------------------------

        // Activate or deactivate measurement noise:
        meas.setNoiseActivated(true);

        // Alternative 1: Measure pellet concentration in one specific cell. Sample point can be freely chosen:
        int[] samplePoint = new int[] {cageDims[0]/2, cageDims[1]/2, cageDims[2]-2};
        measurement = meas.measurePelletConcentrationAtPoint(samplePoint);

        // Alternative 2: Measure total amount of pellets in cage:
        //measurement = meas.measureTotalFeed();

        // Alternative 3: Measure average depth of feed.
        //measurement = meas.measureAverageFeedDepth();


        // -------------------------------------------------------------------------------------------------------------
        // Calculate and return your control signal here:
        // -------------------------------------------------------------------------------------------------------------

        return 50;//nominalFeedRate; // Just keep the initial feeding rate.


    }
}
