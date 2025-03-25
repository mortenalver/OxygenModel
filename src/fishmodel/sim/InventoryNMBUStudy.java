package fishmodel.sim;

import ucar.ma2.ArrayDouble;
import ucar.ma2.InvalidRangeException;
import ucar.nc2.NetcdfFile;
import ucar.nc2.Variable;

import java.io.IOException;
import java.util.Date;

public class InventoryNMBUStudy {

    double[] countVal, weightVal, feedVal;
    public InventoryNMBUStudy(String invDataFile, Date startTime) {
        try {
            NetcdfFile ncfile = NetcdfFile.open(invDataFile);
            Variable time = ncfile.findVariable("time");
            Variable count = ncfile.findVariable("count");
            Variable weight = ncfile.findVariable("weight");
            Variable feed = ncfile.findVariable("feed");

            // Get shape of time variable:
            int[] shape = time.getShape();
            ArrayDouble.D1 tdata = (ArrayDouble.D1) time.read(new int[]{0}, shape);

            Date aTime = null;
            Long[] ltime = new Long[shape[0]];
            int valToUse = -1;
            for (int i=0; i<ltime.length; i++) {
                ltime[i] = (1000*Math.round(tdata.get(i) * 86400 - 7200)); // Subtracting two hours since Java assumes GMT, but it is given in Norwegian summer time.
                aTime = new Date(ltime[i]);
                //System.out.println("i="+i+": "+aTime.toString());
                if (aTime.compareTo(startTime) > 0) {
                    valToUse = Math.max(0, i - 1);
                    break;
                }

            }

            System.out.println("INVENTORY: Date match at index "+valToUse);

            // Get shape of count variable:
            int[] cShape = count.getShape();
            countVal = new double[cShape[1]];
            weightVal = new double[cShape[1]];
            feedVal = new double[cShape[1]];
            ArrayDouble.D2 cdata = (ArrayDouble.D2) count.read(
                    new int[]{valToUse, 0}, new int[] {1, cShape[1]});
            ArrayDouble.D2 wdata = (ArrayDouble.D2) weight.read(
                    new int[]{valToUse, 0}, new int[] {1, cShape[1]});
            ArrayDouble.D2 fdata = (ArrayDouble.D2) feed.read(
                    new int[]{valToUse, 0}, new int[] {1, cShape[1]});
            for (int i=0; i<countVal.length; i++) {
                countVal[i] = cdata.get(0, i);
                weightVal[i] = 1000.*wdata.get(0, i);
                /*if (countVal[i] > 0)
                    feedVal[i] = countVal[i]*weightVal[i]*0.02;*/
                feedVal[i] = 1000.*fdata.get(0, i);
                //System.out.println("Count: "+countVal[i]+", Weight: "+weightVal[i]);
            }

            ncfile.close();
        } catch (IOException e) {
            throw new RuntimeException(e);
        } catch (InvalidRangeException e) {
            throw new RuntimeException(e);
        }
    }
    public  double[] getWeight() {
        return weightVal;
    }

    public double[] getCount() {
        return countVal;
    }

    public double[] getFeed() {
        return feedVal;
    }

}
