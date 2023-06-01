package save;

import ucar.ma2.Array;
import ucar.ma2.ArrayDouble;
import ucar.ma2.DataType;
import ucar.ma2.InvalidRangeException;
import ucar.nc2.Attribute;
import ucar.nc2.Dimension;
import ucar.nc2.NetcdfFileWriteable;
import ucar.nc2.Variable;

import java.io.IOException;
import java.util.ArrayList;

public class SaveForEnKF {

    private static void createTimeDimension(NetcdfFileWriteable ncfile) {
        Dimension tDim = ncfile.addUnlimitedDimension("time");
    }

    private static void createTimeVariable(NetcdfFileWriteable ncfile, String units) {
        ArrayList<Dimension> dimensions = new ArrayList<Dimension>();
        dimensions.add(ncfile.findDimension("time"));
        Variable var = ncfile.addVariable("time", DataType.DOUBLE, dimensions);
        if (units != null)
            var.addAttribute(new Attribute("units", units));
        else
            var.addAttribute(new Attribute("units", "seconds since 2000-01-01 00:00:00"));
    }

    public static NetcdfFileWriteable initializeEnKFFile(String f, int nStates, int N, int m, String timeUnits) {
        try {
            NetcdfFileWriteable ncfile = NetcdfFileWriteable.createNew(f);
            createTimeDimension(ncfile);
            ncfile.setRedefineMode(true);
            ncfile.addDimension("xc", nStates);
            ncfile.addDimension("yc", N);
            ncfile.addDimension("zc", m);

            ncfile.create();
            ncfile.setRedefineMode(true);
            createTimeVariable(ncfile, timeUnits);
            createEnKFVariable(ncfile, "X", "time", "yc", "xc");
            createEnKFVariable(ncfile, "X_a", "time", "yc", "xc");
            createEnKFVariable(ncfile, "X_twin", "time", "xc");
            createEnKFVariable(ncfile, "M", "time", "xc", "zc");
            ncfile.setRedefineMode(false);

            return ncfile;
        } catch (IOException ex) {
            ex.printStackTrace();
            return null;
        }
    }


    public static NetcdfFileWriteable openFile(String f) {
        try {
            NetcdfFileWriteable ncfile = NetcdfFileWriteable.openExisting(f);

            return ncfile;
        } catch (IOException ex) {
            ex.printStackTrace();
            return null;
        }
    }
    public static void createEnKFVariable(NetcdfFileWriteable ncfile, String name, String... dims) {
        try {
            ncfile.setRedefineMode(true);
            ArrayList<Dimension> dimensions = new ArrayList<Dimension>();
            for (String dim : dims) {
                dimensions.add(ncfile.findDimension(dim));
            }
            Variable var = ncfile.addVariable(name, DataType.DOUBLE, dimensions);
            var.addAttribute(new Attribute("scale_factor", 1.0));
            var.addAttribute(new Attribute("add_offset", 0.0));
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public static void save1DVariable(NetcdfFileWriteable ncfile, String name, int tIndex, double[] value, Dimension dim) throws InvalidRangeException, IOException {

        ArrayDouble.D2 varData = new ArrayDouble.D2(1, dim.getLength());

        for (int i1=0; i1<dim.getLength(); i1++)
            varData.set(0, i1, value[i1]);

        ncfile.write(name,  new int[] {tIndex, 0}, varData);
    }

    public static void save2DVariable(NetcdfFileWriteable ncfile, String name, int tIndex, double[][] value, Dimension... dims) throws InvalidRangeException, IOException {

        ArrayDouble.D3 varData = new ArrayDouble.D3(1, dims[0].getLength(), dims[1].getLength());

        for (int i1=0; i1<dims[0].getLength(); i1++)
            for (int i2=0; i2<dims[1].getLength(); i2++)
                varData.set(0, i1, i2, value[i2][i1]);

        ncfile.write(name,  new int[] {tIndex, 0, 0}, varData);
    }

    public static void saveEnKFVariables(NetcdfFileWriteable ncfile, double time, double[][] X, double[][] X_a,
                                         double[] X_twin, double[][] M) {
        try {
            Dimension tDim = ncfile.findDimension("time"),
                    xDim = ncfile.findDimension("xc"),
                    yDim = ncfile.findDimension("yc"),
                    zDim = ncfile.findDimension("zc");
            Array timeData = Array.factory(DataType.DOUBLE, new int[]{1});
            timeData.setDouble(0, time);
            int[] timeOrigin = new int[]{tDim.getLength()};
            ncfile.write("time", timeOrigin, timeData);

            int[] dataOrigin = new int[] {tDim.getLength()-1, 0, 0, 0};

            save2DVariable(ncfile,  "X", tDim.getLength()-1, X, yDim, xDim);
            save2DVariable(ncfile,  "X_a", tDim.getLength()-1, X_a, yDim, xDim);
            save1DVariable(ncfile,  "X_twin", tDim.getLength()-1, X_twin, xDim);
            save2DVariable(ncfile,  "M", tDim.getLength()-1, M, xDim, zDim);
        } catch (InvalidRangeException e) {
            throw new RuntimeException(e);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
}
