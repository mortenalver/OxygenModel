package fishmodel.sim;

import ucar.ma2.ArrayDouble;
import ucar.ma2.InvalidRangeException;
import ucar.nc2.NetcdfFile;
import ucar.nc2.Variable;

import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Random;

public class InputDataNetcdf {

    private Date[] times;
    private Long[] ltime;
    private double[] currentSpeeds, currentDirs, o2ambVal5, o2ambVal10, o2ambVal15,
            feedingVal, temperatures5, temperatures10, temperatures15;
    private double[][] currentSpeeds2, currentDirs2;
    private double currentSpeedMod, currentDirMod;
    private double[] currentSpeedMod2=null, currentDirMod2=null, currentDepths=null;

    long sTime = 0l;
    private int piv = -1;

    private static boolean ADD_CURRENT_NOISE = false;
    private static double CURRENT_COMPONENTS_STD = 0.02;
    private static double currentDirOffset = 25;
    private Random random = new Random();

    private double[] addTemperatureOffsets = new double[] {0, 0, 0};

    private boolean multipleLayers = false; // Set to true if the NetCDF file gives full current profiles
    public InputDataNetcdf(String filepath, boolean useInstantaneousAmbientVals) {

        try {
            NetcdfFile ncfile = NetcdfFile.open(filepath);
            Variable time = ncfile.findVariable("time");
            Variable spd = ncfile.findVariable("extCurrentSpeed");
            Variable direction = ncfile.findVariable("extCurrentDir");
            Variable o2amb5 = ncfile.findVariable(useInstantaneousAmbientVals ? "O2ambient_5" : "O2constAmbient_5");
            Variable o2amb10 = ncfile.findVariable(useInstantaneousAmbientVals ? "O2ambient_10" : "O2constAmbient_10");
            Variable o2amb15 = ncfile.findVariable(useInstantaneousAmbientVals ? "O2ambient_15" : "O2constAmbient_15");
            //Variable feedingBitfield = ncfile.findVariable("feedingBitfield");
            Variable temp5 = ncfile.findVariable("temperature_5");
            Variable temp10 = ncfile.findVariable("temperature_10");
            Variable temp15 = ncfile.findVariable("temperature_15");

            // Get shape of time variable:
            int[] shape = time.getShape();

            // Check if we have 1 series of current data, or multiple layers:
            int[] currentShape = spd.getShape();
            multipleLayers = currentShape.length > 1;

            ArrayDouble.D1 tdata = (ArrayDouble.D1) time.read(new int[]{0}, shape);

            ArrayDouble.D1 csdata=null;
            ArrayDouble.D1 cddata=null;
            ArrayDouble.D2 csdata2=null;
            ArrayDouble.D2 cddata2=null;
            if (multipleLayers) {
                csdata2 = (ArrayDouble.D2) spd.read(new int[]{0, 0}, currentShape);
                cddata2 = (ArrayDouble.D2) direction.read(new int[]{0, 0}, currentShape);
                // Read depth layers for current profiles:
                Variable zc = ncfile.findVariable("zc");
                ArrayDouble.D1 zcAD = (ArrayDouble.D1) zc.read(new int[]{0}, zc.getShape());
                currentDepths = new double[zc.getShape(0)];
                for (int i=0; i<currentDepths.length; i++) {
                    currentDepths[i] = zcAD.get(i);
                }
            } else {
                csdata = (ArrayDouble.D1) spd.read(new int[]{0}, shape);
                cddata = (ArrayDouble.D1) direction.read(new int[]{0}, shape);
            }
            ArrayDouble.D1 o2ambdata5 = (ArrayDouble.D1) o2amb5.read(new int[]{0}, shape);
            ArrayDouble.D1 o2ambdata10 = (ArrayDouble.D1) o2amb10.read(new int[]{0}, shape);
            ArrayDouble.D1 o2ambdata15 = (ArrayDouble.D1) o2amb15.read(new int[]{0}, shape);
            //ArrayDouble.D1 feedingdata = (ArrayDouble.D1) feedingBitfield.read(new int[]{0}, shape);
            ArrayDouble.D1 tempdata5 = (ArrayDouble.D1) temp5.read(new int[]{0}, shape);
            ArrayDouble.D1 tempdata10 = (ArrayDouble.D1) temp10.read(new int[]{0}, shape);
            ArrayDouble.D1 tempdata15 = (ArrayDouble.D1) temp15.read(new int[]{0}, shape);
            times = new Date[shape[0]];
            ltime = new Long[shape[0]];
            if (multipleLayers) {
                currentSpeeds2 = new double[currentShape[0]][currentShape[1]];
                currentDirs2 = new double[currentShape[0]][currentShape[1]];
            } else {
                currentSpeeds = new double[shape[0]];
                currentDirs = new double[shape[0]];
            }
            o2ambVal5 = new double[shape[0]];
            o2ambVal10 = new double[shape[0]];
            o2ambVal15 = new double[shape[0]];
            feedingVal = new double[shape[0]];
            temperatures5 = new double[shape[0]];
            temperatures10 = new double[shape[0]];
            temperatures15 = new double[shape[0]];
            for (int i=0; i<times.length; i++) {
                ltime[i] = (1000*Math.round(tdata.get(i) * 86400 - 7200)); // Subtracting two hours since Java assumes GMT, but it is given in Norwegian summer time.
                times[i] = new Date(ltime[i]);
                if (!multipleLayers) {
                    currentSpeeds[i] = csdata.get(i);
                    currentDirs[i] = cddata.get(i);
                }
                o2ambVal5[i] = o2ambdata5.get(i);
                o2ambVal10[i] = o2ambdata10.get(i);
                o2ambVal15[i] = o2ambdata15.get(i);
                //feedingVal[i] = feedingdata.get(i);
                temperatures5[i] = tempdata5.get(i) + addTemperatureOffsets[0];
                temperatures10[i] = tempdata10.get(i) + addTemperatureOffsets[1];
                temperatures15[i] = tempdata15.get(i) + addTemperatureOffsets[2];
            }
            if (multipleLayers) {
                for (int i=0; i<times.length; i++) {
                    for (int j=0; j<currentShape[1]; j++) {

                        currentSpeeds2[i][j] = csdata2.get(i,j);
                        currentDirs2[i][j] = cddata2.get(i,j);
                    }
                }
            }

            System.out.println("Init file initial time: "+times[0].toString());

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
        System.out.println(startTime.toString());
        System.out.println();
        while (((piv++) < (ltime.length-1)) && (ltime[piv+1]<sTime));
        System.out.println("piv = "+piv);

        computeCurrentWithNoise();
        //System.out.println(ltime[piv]);
    }

    public boolean advance(double time) {
        int oldPiv = piv;
        long newTime = sTime + (long)(time*1000);
        while ((piv < ltime.length-2) && ltime[piv+1]<newTime)
            piv++;
        if (piv == ltime.length-2)
            System.out.println("No more observation data available!");

        if (piv > oldPiv)
            computeCurrentWithNoise();

        return (piv > oldPiv);
    }

    private void computeCurrentWithNoise() {
        if (multipleLayers) {
            computeCurrentProfileWithNoise();
            return;
        }
        if (ADD_CURRENT_NOISE) {
            //System.out.println("Before: speed = "+currentSpeeds[piv]+", dir = "+currentDirs[piv]);
            // Pick normally distributed random numbers for x and y components:
            double xnoise = random.nextGaussian()*CURRENT_COMPONENTS_STD,
                    ynoise = random.nextGaussian()*CURRENT_COMPONENTS_STD;
            // Decompose measured current vector and add noise values:
            double xcomp = currentSpeeds[piv]*Math.cos(currentDirs[piv]*Math.PI/180.),
                    ycomp = currentSpeeds[piv]*Math.sin(currentDirs[piv]*Math.PI/180.);
            //System.out.println("xcomp = "+xcomp+", ycomp = "+ycomp);
            double xcompMod =  xcomp + xnoise,
                    ycompMod = ycomp + ynoise;
            //System.out.println("xcompmod = "+xcompMod+", ycompmod = "+ycompMod);
            // Compute speed and direction with noise:
            currentSpeedMod = Math.sqrt(xcompMod*xcompMod + ycompMod*ycompMod);
            currentDirMod = (180./Math.PI)*Math.atan2(ycompMod, xcompMod);
            if (currentDirMod < 0)
                currentDirMod += 360.;
            //System.out.println("After: speed = "+currentSpeedMod+", dir  "+currentDirMod);
        }
        else {
            currentSpeedMod = currentSpeeds[piv];
            currentDirMod = currentDirs[piv] + currentDirOffset;
        }
    }

    private void computeCurrentProfileWithNoise() {
        if (currentSpeedMod2 == null) {
            System.out.println("Initializing profile: "+currentSpeeds2[0].length);
            currentSpeedMod2 = new double[currentSpeeds2[0].length];
            currentDirMod2 = new double[currentSpeeds2[0].length];
        }
        for (int i=0; i<currentSpeeds2[0].length; i++) {
            currentSpeedMod2[i] = currentSpeeds2[piv][i];
            currentDirMod2[i] = currentDirs2[piv][i];
        }
    }

    public double getTemperature5() { return temperatures5[piv]; }
    public double getTemperature10() { return temperatures10[piv]; }
    public double getTemperature15() { return temperatures15[piv]; }

    public double getExtCurrentSpeed() {
        return currentSpeedMod; //currentSpeeds[piv];
    }

    public double[] getExtCurrentSpeedProfile() {
        return currentSpeedMod2;
    }

    public double getExtCurrentDir() {
        return currentDirMod;//currentDirs[piv];
    }

    public double[] getExtCurrentDirProfile() {
        return currentDirMod2;
    }

    public double[] getCurrentDepths() {
        return currentDepths;
    }

    public double getO2Ambient5() {
        return o2ambVal5[piv];
    }
    public double getO2Ambient10() {
        return o2ambVal10[piv];
    }
    public double getO2Ambient15() {
        return o2ambVal15[piv];
    }

    public double getFeedingRate() {
        return feedingVal[piv];
    }
}
