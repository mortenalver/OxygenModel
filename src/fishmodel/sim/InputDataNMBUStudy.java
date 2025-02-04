package fishmodel.sim;

import ucar.ma2.ArrayDouble;
import ucar.ma2.InvalidRangeException;
import ucar.nc2.NetcdfFile;
import ucar.nc2.Variable;

import java.io.IOException;
import java.util.Date;
import java.util.Random;

public class InputDataNMBUStudy {

    private Date[] times;
    private Long[] ltime;

    private double[][] currentSpeeds, currentDirs, temps;
    private double[] currentSpeedMod =null, currentDirMod =null, tempMod=null, currentDepths=null;

    long sTime = 0l;
    private int piv = -1;

    private static boolean ADD_CURRENT_NOISE = false;
    private static double CURRENT_COMPONENTS_STD = 0.02;
    private static double currentDirOffset = 25;
    private Random random = new Random();

    private double[] addTemperatureOffsets = new double[] {0, 0, 0};


    public InputDataNMBUStudy(String filepath) {

        try {
            System.out.println("Opening input file: "+filepath);
            NetcdfFile ncfile = NetcdfFile.open(filepath);
            Variable time = ncfile.findVariable("time");
            Variable spd = ncfile.findVariable("extCurrentSpeed");
            Variable direction = ncfile.findVariable("extCurrentDir");
            //Variable o2amb5 = ncfile.findVariable(useInstantaneousAmbientVals ? "O2ambient_5" : "O2constAmbient_5");
            Variable temp = ncfile.findVariable("temperature");

            // Get shape of time variable:
            int[] shape = time.getShape();

            // Check if we have 1 series of current data, or multiple layers:
            int[] currentShape = spd.getShape();

            ArrayDouble.D1 tdata = (ArrayDouble.D1) time.read(new int[]{0}, shape);
            ArrayDouble.D2 csdata = (ArrayDouble.D2) spd.read(new int[]{0, 0}, currentShape);
            ArrayDouble.D2 cddata = (ArrayDouble.D2) direction.read(new int[]{0, 0}, currentShape);
            ArrayDouble.D2 tempdata = (ArrayDouble.D2) temp.read(new int[]{0, 0}, currentShape);

            // Read depth layers for current profiles:
            Variable zc = ncfile.findVariable("zc");
            ArrayDouble.D1 zcAD = (ArrayDouble.D1) zc.read(new int[]{0}, zc.getShape());
            currentDepths = new double[zc.getShape(0)];
            for (int i=0; i<currentDepths.length; i++) {
                currentDepths[i] = zcAD.get(i);
            }

            times = new Date[shape[0]];
            ltime = new Long[shape[0]];
            currentSpeeds = new double[currentShape[0]][currentShape[1]];
            currentDirs = new double[currentShape[0]][currentShape[1]];
            temps = new double[currentShape[0]][currentShape[1]];
            for (int i=0; i<times.length; i++) {
                ltime[i] = (1000*Math.round(tdata.get(i) * 86400 - 7200)); // Subtracting two hours since Java assumes GMT, but it is given in Norwegian summer time.
                times[i] = new Date(ltime[i]);

            }
            for (int i=0; i<times.length; i++) {
                for (int j=0; j<currentShape[1]; j++) {

                    currentSpeeds[i][j] = csdata.get(i,j);
                    currentDirs[i][j] = cddata.get(i,j);
                    temps[i][j] = tempdata.get(i,j);
                }
            }


            ncfile.close();
        } catch (
                IOException ex) {
            ex.printStackTrace();

        } catch (
                InvalidRangeException e) {
            e.printStackTrace();

        }
    }

    public void setStartTime(Date startTime) {
        sTime = startTime.getTime();
        System.out.println("StartTime: "+startTime.toString());
        System.out.println();
        while (((piv++) < (ltime.length-1)) && (ltime[piv+1]<sTime));

        computeCurrent();

    }

    public boolean advance(double time) {
        int oldPiv = piv;
        long newTime = sTime + (long)(time*1000);
        while ((piv < ltime.length-2) && ltime[piv+1]<newTime)
            piv++;
        if (piv == ltime.length-2)
            System.out.println("No more observation data available!");

        if (piv > oldPiv)
            computeCurrent();

        return (piv > oldPiv);
    }



    private void computeCurrent() {
        if (currentSpeedMod == null) {
            currentSpeedMod = new double[currentSpeeds[0].length];
            currentDirMod = new double[currentSpeeds[0].length];
            tempMod = new double[currentSpeeds[0].length];
        }
        for (int i = 0; i< currentSpeeds[0].length; i++) {
            currentSpeedMod[i] = currentSpeeds[piv][i];
            currentDirMod[i] = currentDirs[piv][i];
            //System.out.println("i="+i+", speed="+currentSpeedMod[i]+", dir="+currentDirMod[i]);
            tempMod[i] = temps[piv][i];
        }
    }

    public double[] getTemperatureProfile() { return tempMod; }


    public double[] getExtCurrentSpeedProfile() {
        return currentSpeedMod;
    }

    public double[] getExtCurrentDirProfile() {
        return currentDirMod;
    }

    public double[] getCurrentDepths() {
        return currentDepths;
    }

    public double getO2Ambient5() {
        return 9.;
    }

    public double getO2Ambient10() {
        return 9.;
    }

    public double getO2Ambient15() {
        return 9.;
    }
}
