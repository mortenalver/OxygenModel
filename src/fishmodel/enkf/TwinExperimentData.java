package fishmodel.enkf;

import ucar.ma2.ArrayDouble;
import ucar.ma2.InvalidRangeException;
import ucar.nc2.NetcdfFile;
import ucar.nc2.Variable;

import java.io.IOException;

/**
 * This class is used to read O2 sensor values from a previous output file in order to provide measurements for
 * twin experiment.
 */
public class TwinExperimentData {

    private String[] varNames = new String[] {"C_5", "C_10", "C_15",
            "M1_5", "M1_10", "M1_15",
            "M2_5", "M2_10", "M2_15",
            "M3_5", "M3_10", "M3_15"
    };

    private double[] time;
    private double[][] o2;

    int piv = 0; // Internal value for the index of last returned value.


    public TwinExperimentData(String filepath) {
        try {
            NetcdfFile ncfile = NetcdfFile.open(filepath);
            Variable vtime = ncfile.findVariable("time");
            int[] shape = vtime.getShape();
            for (int i = 0; i < shape.length; i++) {
                int i1 = shape[i];
                System.out.println(i1);
            }
            time = new double[shape[0]];
            ArrayDouble.D1 tVal = (ArrayDouble.D1)vtime.read(new int[] {0}, shape);
            for (int i=0; i<time.length; i++)
                time[i] = tVal.get(i);

            o2 = new double[varNames.length][time.length];
            for (int i=0; i<varNames.length; i++) {
                Variable o2Var = ncfile.findVariable(varNames[i]);
                shape = o2Var.getShape();
                ArrayDouble.D1 o2Val = (ArrayDouble.D1)o2Var.read(new int[] {0}, shape);
                for (int j=0; j<time.length; j++) {
                    o2[i][j] = o2Val.get(j);
                }
            }

        } catch (IOException e) {
            throw new RuntimeException(e);
        } catch (InvalidRangeException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Get the first measurement that has time stamp equal to or smaller than t. The value of t is assumed to
     * increase monotonically over repeated calls - so we never need to go back to earlier measurement times.
     * Times are given in seconds since the start of a simulation, so absolute time values are not considered.
     * @param t The time for which we should find the measurement.
     * @return An array of measurements, or null if not found.
     */
    public double[] getMeasurements(double t) {
        // Count upwards until we have the correct time:
        while ((time[piv] < t) && (piv < time.length-1))
            piv++;
        if (piv < time.length) {
            double[] res = new double[varNames.length];
            for (int i=0; i<varNames.length; i++)
                res[i] = o2[i][piv];
            return res;
        }
        else {
            System.out.println("TwinExperimentData: measurement time not found.");
            return null;
        }
    }
}
