package fishmodel.sim;

/**
 * Created by malv on 12.06.15.
 */
public class StoredState {

    private SimulationInputs si;
    private double[][][] fc;

    /**
     * Constructs a stored state instance for the given simulation inputs and model state.
     * The objects are copied, not stored by reference.
     * @param inputs the simulation inputs
     * @param feed the feed concentration
     */
    public StoredState(SimulationInputs inputs, double[][][] feed) {
        this.si = (SimulationInputs)inputs.clone();
        fc = new double[feed.length][feed[0].length][feed[0][0].length];
        copyArray(feed, fc);
    }

    public void restoreState(double[][][] feed) {
        copyArray(fc, feed);
    }

    protected void copyArray(double[][][] from, double[][][] to) {
        for (int i=0; i<from.length; i++) {
            for (int j=0; j<from[i].length; j++) {
                for (int k=0; k<from[i][j].length; k++)
                    to[i][j][k] = from[i][j][k];
            }
        }
    }

    public SimulationInputs getSimulationInputs() {
        return (SimulationInputs)si.clone();
    }

    public double[][][] getFeedConcentration() {
        return fc;
    }
}
