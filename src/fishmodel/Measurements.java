package fishmodel;

import fishmodel.enkf.Util;
import org.ejml.data.DMatrixRMaj;

import java.util.HashMap;
import java.util.LinkedHashMap;

public class Measurements {

    public static class MeasurementSet {
        public String[] names;
        public int[][] pos;

        public double std = 0;
    }

    public static MeasurementSet setupSensorPositionsBjoroya(int[] cageDims, double dxy, double dz, double rad) {
        MeasurementSet ms = new MeasurementSet();
        ms.std = 0.025;

        ms.names = new String[] {"C_5", "C_10", "C_15", "M1_5", "M1_10", "M1_15", "M2_5", "M2_10", "M2_15",
                "M3_5", "M3_10", "M3_15"};
                //"oM1_5", "oM1_10", "oM1_15", "oM2_5", "oM2_10", "oM2_15", "oM3_5", "oM3_10", "oM3_15"};

        // angles: M1 128.2948, M2 2.8445, M3 246.8427
        double angle1 = 128.2948, angle2 = 2.8445, angle3 = 246.8427;
        double[] sensorAngles = new double[] {0, 0, 0, angle1, angle1, angle1, angle2, angle2, angle2,
                angle3, angle3, angle3,
                angle1, angle1, angle1, angle2, angle2, angle2,
                angle3, angle3, angle3};
        //double[] o2Rad = new double[] {0.6, 0.6, 0.6, 0.6}; // Distance from centre as fraction of cage radius
        //double[] o2Rad = new double[] {0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1}; // Distance from centre as fraction of cage radius
        double[] o2Rad = new double[] {0, 0, 0, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.95, 0.95, 0.9};
                //1.02, 1.02, 1.02, 1.05, 1.05, 1.05, 1.01, 1.01, 1.01}; // Distance from centre as fraction of cage radius
        double[] o2Depth = new double[] {5, 10, 15, 5, 10, 15, 5, 10, 15, 5, 10, 15};
                //5, 10, 15, 5, 10, 15, 5, 10, 15}; // Sensor depth (m)

        // Set up o2 sensor positions in grid:
        ms.pos = new int[ms.names.length][3];
        for (int i=0; i<ms.names.length; i++) {
            int xDist = (int)(o2Rad[i]*rad*Math.sin(Math.PI*sensorAngles[i]/180.)/dxy + cageDims[0]/2);
            int yDist = (int)(o2Rad[i]*rad*Math.cos(Math.PI*sensorAngles[i]/180.)/dxy + cageDims[1]/2);
            int zDist = (int)(o2Depth[i]/dz);
            ms.pos[i][0] = xDist;
            ms.pos[i][1] = yDist;
            ms.pos[i][2] = zDist;
            //System.out.println(ms.names[i]+": "+ms.pos[i][0]+" , "+ms.pos[i][1]+" , "+ms.pos[i][2]);
        }

        return ms;
    }


    public static DMatrixRMaj getMeasurementModel(int[] cageDims, MeasurementSet ms) {
        int n = cageDims[0]*cageDims[1]*cageDims[2];
        DMatrixRMaj M = new DMatrixRMaj(ms.names.length, n);
        for (int i=0; i<ms.names.length; i++) {
            int stateNum = Util.getStateIndex(ms.pos[i][0], ms.pos[i][1], ms.pos[i][2], cageDims);
            //System.out.println("meas "+i+": state "+stateNum);
            M.set(i, stateNum, 1);
        }
        return M;
    }
}
