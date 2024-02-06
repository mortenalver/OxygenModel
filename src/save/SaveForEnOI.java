package save;

import fishmodel.Measurements;
import ucar.ma2.*;
import ucar.nc2.Attribute;
import ucar.nc2.Dimension;
import ucar.nc2.NetcdfFileWriteable;
import ucar.nc2.Variable;

import java.io.IOException;
import java.util.ArrayList;

public class SaveForEnOI {

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

    public static NetcdfFileWriteable initializeFile(String f, int nStates, int N, int m, int m_all, int nPar, String timeUnits,
                                                     int[] cageDims, double dxy, Measurements.MeasurementSet ms) {
        try {
            NetcdfFileWriteable ncfile = NetcdfFileWriteable.createNew(f);
            addMeasurementAttributes(ncfile, ms);
            ncfile.addGlobalAttribute("nParameters", nPar);
            createTimeDimension(ncfile);
            ncfile.setRedefineMode(true);
            ncfile.addDimension("xc", nStates);
            ncfile.addDimension("yc", N);
            ncfile.addDimension("zc", m);
            ncfile.addDimension("zc2", m_all);

            ncfile.create();
            ncfile.setRedefineMode(true);
            createTimeVariable(ncfile, timeUnits);
            createVariable(ncfile, "X_static", "yc", "xc");
            createVariable(ncfile, "x_f", "time", "xc");
            createVariable(ncfile, "x_a", "time", "xc");
            createVariable(ncfile, "M", "time", "xc", "zc");
            createVariable(ncfile, "deviation", "time", "zc");
            createVariable(ncfile, "deviation_all", "time", "zc2");
            createVariable(ncfile, "K", "zc", "xc");
            createVariable(ncfile, "Xloc", "zc", "xc");

            ArrayInt cdA = new ArrayInt(new int[] {cageDims.length});
            for (int i = 0; i < cageDims.length; i++) {
                cdA.setInt(i, cageDims[i]);
            }

            ncfile.addGlobalAttribute("cageDims", cdA);
            ncfile.addGlobalAttribute("dxy", dxy);
            ncfile.setRedefineMode(false);

            return ncfile;
        } catch (IOException ex) {
            ex.printStackTrace();
            return null;
        }
    }

    private static void addMeasurementAttributes(NetcdfFileWriteable ncfile, Measurements.MeasurementSet ms) {
        try {
            ncfile.setRedefineMode(true);
            // Store delimited string with sensor names:
            String[] names = ms.names;
            StringBuilder sb = new StringBuilder();
            for (int i=0; i<names.length; i++) {
                sb.append(names[i]);
                if (i<names.length-1)
                    sb.append("\t");
            }
            ncfile.addGlobalAttribute("sensorNames", sb.toString());
            // Store sensor positions (x, y and z):
            String chrs = "XYZ";
            for (int i=0; i<3; i++) {
                ArrayInt.D1 positions = new ArrayInt.D1(ms.pos.length);
                for (int j=0; j<ms.pos.length; j++)
                    positions.set(j, ms.pos[j][i]);
                ncfile.addGlobalAttribute("sensorPositions"+chrs.substring(i,i+1), positions);
            }
            // Store which sensors are assimilated:
            int[] as = Measurements.getSensorsToAssimilate();
            ArrayInt.D1 assimSensors = new ArrayInt.D1(as.length);
            for (int i = 0; i < as.length; i++) {
                assimSensors.set(i, as[i]);
            }
            ncfile.addGlobalAttribute("sensorsAssimilated", assimSensors);

        } catch (IOException e) {
            e.printStackTrace();
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
    public static void createVariable(NetcdfFileWriteable ncfile, String name, String... dims) {
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

    public static void save2DVariableNoTime(NetcdfFileWriteable ncfile, String name, double[][] value, Dimension... dims) throws InvalidRangeException, IOException {

        ArrayDouble.D2 varData = new ArrayDouble.D2(dims[0].getLength(), dims[1].getLength());

        for (int i1=0; i1<dims[0].getLength(); i1++)
            for (int i2=0; i2<dims[1].getLength(); i2++)
                varData.set(i1, i2, value[i2][i1]);

        ncfile.write(name,  new int[] {0, 0}, varData);
    }

    public static void save2DVariable(NetcdfFileWriteable ncfile, String name, int tIndex, double[][] value, Dimension... dims) throws InvalidRangeException, IOException {

        ArrayDouble.D3 varData = new ArrayDouble.D3(1, dims[0].getLength(), dims[1].getLength());

        for (int i1=0; i1<dims[0].getLength(); i1++)
            for (int i2=0; i2<dims[1].getLength(); i2++)
                varData.set(0, i1, i2, value[i2][i1]);

        ncfile.write(name,  new int[] {tIndex, 0, 0}, varData);
    }

    public static void saveStaticVariables(NetcdfFileWriteable ncfile, double time, double[][] X_static,
                                     double[][] Xloc, double[][] K) {
        try {
            Dimension xDim = ncfile.findDimension("xc"),
                    yDim = ncfile.findDimension("yc"),
                    zDim = ncfile.findDimension("zc"),
                    zDim2 = ncfile.findDimension("zc2");

            save2DVariableNoTime(ncfile,  "X_static", X_static, yDim, xDim);
            save2DVariableNoTime(ncfile,  "K", K, zDim, xDim);
            save2DVariableNoTime(ncfile,  "Xloc", Xloc, zDim, xDim);
        } catch (InvalidRangeException e) {
            throw new RuntimeException(e);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }


    public static void saveVariables(NetcdfFileWriteable ncfile, double time,
                                     double[] x_f, double[] x_a, double[] dev,
                                     double[] dev_all) {
        try {
            Dimension tDim = ncfile.findDimension("time"),
                    xDim = ncfile.findDimension("xc"),
                    yDim = ncfile.findDimension("yc"),
                    zDim = ncfile.findDimension("zc"),
                    zDim2 = ncfile.findDimension("zc2");
            Array timeData = Array.factory(DataType.DOUBLE, new int[]{1});
            timeData.setDouble(0, time);
            int[] timeOrigin = new int[]{tDim.getLength()};
            ncfile.write("time", timeOrigin, timeData);


            //save2DVariable(ncfile,  "X_a", tDim.getLength()-1, X_a, yDim, xDim);

            //save2DVariable(ncfile,  "M", tDim.getLength()-1, M, xDim, zDim);
            save1DVariable(ncfile,  "x_f", tDim.getLength()-1, x_f, xDim);
            save1DVariable(ncfile,  "x_a", tDim.getLength()-1, x_a, xDim);
            save1DVariable(ncfile,  "deviation", tDim.getLength()-1, dev, zDim);
            save1DVariable(ncfile,  "deviation_all", tDim.getLength()-1, dev_all, zDim2);
            //System.out.println("K double: "+K.length+" by "+K[0].length);
            //System.out.println("xdim: "+xDim.getLength());
            //System.out.println("zdim: "+zDim.getLength());
            //save2DVariable(ncfile,  "K", tDim.getLength()-1, K, zDim, xDim);
            //save2DVariable(ncfile,  "Xloc", tDim.getLength()-1, Xloc, zDim, xDim);

        } catch (InvalidRangeException e) {
            throw new RuntimeException(e);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
}
