package fishmodel;

import fishmodel.enkf.Util;
import org.ejml.data.DMatrixRMaj;

import java.util.ArrayList;
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
        ms.std = 0.1;

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

    public static MeasurementSet setupSensorPositionsFullFarm(int[] cageDims, double dxy, double dz, double rad,
                                                              double farmRotation, boolean includeExtPos,
                                                              double frameSize,
                                                              ArrayList<double[]> cagePositions, double[] mainCageCenter) {
        MeasurementSet ms = new MeasurementSet();
        ms.std = 0.1;

        ms.names = new String[] {"C_5", "C_10", "C_15", "M1_5", "M1_10", "M1_15", "M2_5", "M2_10", "M2_15",
                "M3_5", "M3_10", "M3_15"};

        int nCageSensors = ms.names.length;

        double[] o2Depth = new double[] {5, 10, 15, 5, 10, 15, 5, 10, 15, 5, 10, 15}; // Sensor depth (m)
        if (includeExtPos) { // Add sensor names for the external sensor position:
            String[] tmp = ms.names;
            double[] tmpD = o2Depth;
            ms.names = new String[tmp.length+3];
            o2Depth = new double[tmpD.length+3];
            for (int i=0; i<tmp.length; i++) {
                ms.names[i] = tmp[i];
                o2Depth[i] = tmpD[i];
            }
            ms.names[tmp.length] = "Ext_5";
            ms.names[tmp.length+1] = "Ext_10";
            ms.names[tmp.length+2] = "Ext_15";
            o2Depth[tmpD.length] = 5;
            o2Depth[tmpD.length+1] = 10;
            o2Depth[tmpD.length+2] = 15;
        }

        // angles: M1 128.2948, M2 2.8445, M3 246.8427
        double angle1 = 128.2948-farmRotation, angle2 = 2.8445-farmRotation, angle3 = 246.8427-farmRotation; // Angles are adjusted for domain rotation
        double[] sensorAngles = new double[] {0, 0, 0, angle1, angle1, angle1, angle2, angle2, angle2,
                angle3, angle3, angle3};
        //double[] o2Rad = new double[] {0.6, 0.6, 0.6, 0.6}; // Distance from centre as fraction of cage radius
        //double[] o2Rad = new double[] {0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1}; // Distance from centre as fraction of cage radius
        double[] o2Rad = new double[] {0, 0, 0, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.91, 0.91, 0.91};

        // Set up o2 sensor positions in grid:
        ms.pos = new int[ms.names.length][3];
        for (int i=0; i<ms.names.length; i++) {
            if (i < nCageSensors) {
                int xDist = (int) (o2Rad[i] * rad * Math.sin(Math.PI * sensorAngles[i] / 180.) / dxy + mainCageCenter[0] / dxy);
                int yDist = (int) (o2Rad[i] * rad * Math.cos(Math.PI * sensorAngles[i] / 180.) / dxy + mainCageCenter[1] / dxy);
                int zDist = (int) (o2Depth[i] / dz);
                System.out.println(ms.names[i]+": "+xDist+" , "+yDist+" , "+zDist);
                ms.pos[i][0] = xDist;
                ms.pos[i][1] = yDist;
                ms.pos[i][2] = zDist;
            } else {
                // Not one of the cage sensors, so this is the external position:
                // Place relative to the second cage:
                double[] cagePos0 = cagePositions.get(1);
                int xDist = (int)((cagePos0[0] - frameSize*2.39)/dxy);
                int yDist = (int)((cagePos0[1] - frameSize*0.36)/dxy);
                int zDist = (int) (o2Depth[i] / dz);
                System.out.println(ms.names[i] + ": " + xDist + " , " + yDist + " , " + zDist);
                ms.pos[i][0] = xDist;
                ms.pos[i][1] = yDist;
                ms.pos[i][2] = zDist;
            }

        }
    return ms;
    }

    public static int[] getSensorsToAssimilate() {
        return getSensorsToAssimilateBjoroya();
    }

    public static int[] getSensorsToAssimilateBjoroya() {

        //return new int[] {0, 1, 2}; // Centre only
        //return new int[] {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}; // All sensors
        //return new int [] {0, 1, 2, 4, 7, 10};  // All at 10 m and all in centre
        //return new int[] {0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}; // All sensors except C10
        return new int[] {3, 4, 5, 6, 7, 8, 9, 10, 11}; // All M sensors

        // Copied from python code:
        //return (0, 3, 6, 9) # All at 5 m
        //return (1, 4, 7, 10)  # All at 10 m
        //return (4, 7, 10)  # Ring measurements at 10 m


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

    public static DMatrixRMaj getMeasurementModelBjoroya(int[] cageDims, int nPar, MeasurementSet ms, boolean allSensors) {
        int[] activeSensors;
        if (allSensors) { // Return model for all sensors
            activeSensors = new int[ms.names.length];
            for (int i=0; i< activeSensors.length; i++)
                activeSensors[i] = i;
        } else { // Return model for active sensors only
            activeSensors = getSensorsToAssimilateBjoroya();
        }
        int n = cageDims[0]*cageDims[1]*cageDims[2];
        DMatrixRMaj M = new DMatrixRMaj(activeSensors.length, n + nPar);
        for (int i=0; i<activeSensors.length; i++) {
            int stateNum = Util.getStateIndex(ms.pos[activeSensors[i]][0],
                    ms.pos[activeSensors[i]][1], ms.pos[activeSensors[i]][2], cageDims);
            //System.out.println("meas "+i+": state "+stateNum);
            M.set(i, stateNum, 1);
        }
        return M;
    }
}
