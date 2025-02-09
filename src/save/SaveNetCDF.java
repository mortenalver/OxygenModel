package save;

import fishmodel.CageStats;
import fishmodel.Measurements;
import ucar.ma2.*;
import ucar.nc2.Attribute;
import ucar.nc2.Dimension;
import ucar.nc2.NetcdfFileWriteable;
import ucar.nc2.Variable;

import java.io.IOException;
import java.util.ArrayList;

/**
 *
 */
public class SaveNetCDF {

    public static NetcdfFileWriteable initializeFile(String f, int[] cageDims, double dxy, double dz) {
        return initializeFile(f, cageDims, dxy, dz, null, null);
    }


    public static NetcdfFileWriteable initializeFile(String f, int[] cageDims, double dxy, double dz, String timeUnits) {
        return initializeFile(f, cageDims, dxy, dz, timeUnits, null);
    }
    public static NetcdfFileWriteable initializeFile(String f, int[] cageDims, double dxy, double dz, String timeUnits,
                                                     Measurements.MeasurementSet ms) {
        try {
            NetcdfFileWriteable ncfile = NetcdfFileWriteable.createNew(f);
            createTimeDimension(ncfile);
            createCageDimensions(ncfile, cageDims, dxy, dz);
            addMatdispAttributes(ncfile, dxy, dz);
            if (ms != null)
                addMeasurementAttributes(ncfile, ms);
            ncfile.create();
            ncfile.setRedefineMode(true);
            createTimeVariable(ncfile, timeUnits);


            return ncfile;
        } catch (IOException ex) {
            ex.printStackTrace();
            return null;
        }
    }


    public static NetcdfFileWriteable initializeFishDataFile(String f, int fishCount) {
        try {
            NetcdfFileWriteable ncfile = NetcdfFileWriteable.createNew(f);
            createTimeDimension(ncfile);
            ncfile.addDimension("xc", fishCount);
            addMatdispAttributes(ncfile, 1, 1);
            ncfile.create();
            ncfile.setRedefineMode(true);
            createTimeVariable(ncfile, null);
            return ncfile;
        } catch (IOException ex) {
            ex.printStackTrace();
            return null;
        }
    }

    /*public static NetcdfFileWriteable initializeScalarDataFile(String f, int fishCount) {
        try {
            NetcdfFileWriteable ncfile = NetcdfFileWriteable.createNew(f);
            createTimeDimension(ncfile);
            ncfile.addDimension("xc", fishCount);
            addMatdispAttributes(ncfile, 1, 1);
            ncfile.create();
            ncfile.setRedefineMode(true);
            createTimeVariable(ncfile);
            return ncfile;
        } catch (IOException ex) {
            ex.printStackTrace();
            return null;
        }
    }*/

    private static void createTimeDimension(NetcdfFileWriteable ncfile) {
        Dimension tDim = ncfile.addUnlimitedDimension("time");
    }

    private static void createCageDimensions(NetcdfFileWriteable ncfile, int[] dims, double dxy, double dz) {
        try {
            ncfile.setRedefineMode(true);
            ncfile.addDimension("xc", dims[0]);
            ncfile.addDimension("yc", dims[1]);
            ncfile.addDimension("zc", dims[2]);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void setGlobalAttribute(NetcdfFileWriteable ncfile, String attName, double value) {
        try {
            ncfile.setRedefineMode(true);
            ncfile.addGlobalAttribute(attName, value);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void setGlobalAttribute(NetcdfFileWriteable ncfile, String attName, int value) {
        try {
            ncfile.setRedefineMode(true);
            ncfile.addGlobalAttribute(attName, value);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void setGlobalAttribute(NetcdfFileWriteable ncfile, String attName, String value) {
        try {
            ncfile.setRedefineMode(true);
            ncfile.addGlobalAttribute(attName, value);
        } catch (IOException e) {
            e.printStackTrace();
        }
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

    public static void createCageVariables(NetcdfFileWriteable ncfile, String... name) {
        for (int i = 0; i < name.length; i++) {
            String s = name[i];
            createCageVariable(ncfile, s);
        }
    }

    public static void createScalarVariables(NetcdfFileWriteable ncfile, String... name) {
        for (int i = 0; i < name.length; i++) {
            String s = name[i];
            createScalarVariable(ncfile, s);
        }
    }

    public static void createCageVariable(NetcdfFileWriteable ncfile, String name) {
        try {
            ncfile.setRedefineMode(true);
            ArrayList<Dimension> dimensions = new ArrayList<Dimension>();
            dimensions.add(ncfile.findDimension("time"));
            dimensions.add(ncfile.findDimension("zc"));
            dimensions.add(ncfile.findDimension("yc"));
            dimensions.add(ncfile.findDimension("xc"));
            Variable var = ncfile.addVariable(name, DataType.DOUBLE, dimensions);
            var.addAttribute(new Attribute("scale_factor", 1.0));
            var.addAttribute(new Attribute("add_offset", 0.0));
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public static void createProfileVariable(NetcdfFileWriteable ncfile, String name, int dim) {
        try {
            ncfile.setRedefineMode(true);
            ArrayList<Dimension> dimensions = new ArrayList<Dimension>();
            dimensions.add(ncfile.findDimension("time"));
            // Second dimension:
            if (dim==0)
                dimensions.add(ncfile.findDimension("xc"));
            else if (dim==1)
                dimensions.add(ncfile.findDimension("yc"));
            else
                dimensions.add(ncfile.findDimension("zc"));

            Variable var = ncfile.addVariable(name, DataType.DOUBLE, dimensions);
            var.addAttribute(new Attribute("scale_factor", 1.0));
            var.addAttribute(new Attribute("add_offset", 0.0));
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public static void createFishDataVariable(NetcdfFileWriteable ncfile, String name) {
        try {
            ncfile.setRedefineMode(true);
            ArrayList<Dimension> dimensions = new ArrayList<Dimension>();
            dimensions.add(ncfile.findDimension("time"));
            dimensions.add(ncfile.findDimension("xc"));
            Variable var = ncfile.addVariable(name, DataType.DOUBLE, dimensions);
            var.addAttribute(new Attribute("scale_factor", 1.0));
            var.addAttribute(new Attribute("add_offset", 0.0));
        } catch (IOException e) {
            e.printStackTrace();
        }

    }
    public static void createScalarVariable(NetcdfFileWriteable ncfile, String name) {
        try {
            ncfile.setRedefineMode(true);
            ArrayList<Dimension> dimensions = new ArrayList<Dimension>();
            dimensions.add(ncfile.findDimension("time"));
            Variable var = ncfile.addVariable(name, DataType.DOUBLE, dimensions);
            var.addAttribute(new Attribute("scale_factor", 1.0));
            var.addAttribute(new Attribute("add_offset", 0.0));
        } catch (IOException e) {
            e.printStackTrace();
        }

    }


    private static void addMatdispAttributes(NetcdfFileWriteable ncfile, double dxy, double dz) {
        try {
            ncfile.setRedefineMode(true);
            ncfile.addGlobalAttribute("TITLE", "AQUAEXCEL");
            ncfile.addGlobalAttribute("_FillValue", -32768);
            ncfile.addGlobalAttribute("grid_mapping_name", "polar_stereographic");
            ncfile.addGlobalAttribute("horizontal_resolution", dxy);
            ncfile.addGlobalAttribute("vertical_resolution", dz);
            ArrayDouble.D1 northPole = new ArrayDouble.D1(2);
            northPole.set(0, 0);
            northPole.set(1, 0);
            ncfile.addGlobalAttribute("coordinate_north_poke", northPole);
            ncfile.addGlobalAttribute("latitude_of_projection_origin", 90);
            ncfile.addGlobalAttribute("standard_parallel", 0.0);
        } catch (IOException e) {
            e.printStackTrace();
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

    public static void saveCageVariable(NetcdfFileWriteable ncfile, double time, String name, double[][][] value, boolean[][][] mask, boolean writeTime) {
        ArrayList<String> names = new ArrayList<String>();
        names.add(name);
        ArrayList<double[][][]> values = new ArrayList<double[][][]>();
        values.add(value);
        ArrayList<boolean[][][]> masks = null;
        if (mask != null) {
            masks = new ArrayList<boolean[][][]>();
            masks.add(mask);
        }
        saveCageVariables(ncfile, time, names, values, masks, writeTime);
    }

    public static void saveCageVariables(NetcdfFileWriteable ncfile, double time, ArrayList<String> names, ArrayList<double[][][]> values,
                                         ArrayList<boolean[][][]> masks, boolean writeTime) {
        try {
            ncfile.setRedefineMode(false);
            Dimension tDim = ncfile.findDimension("time"),
                    xDim = ncfile.findDimension("xc"),
                    yDim = ncfile.findDimension("yc"),
                    zDim = ncfile.findDimension("zc");

            if (writeTime) {
                Array timeData = Array.factory(DataType.DOUBLE, new int[] {1});
                timeData.setDouble(0, time);
                int[] timeOrigin = new int[] {tDim.getLength()};
                ncfile.write("time", timeOrigin, timeData);
            }
            int[] dataOrigin = new int[] {tDim.getLength()-1, 0, 0, 0};

            for (int i=0; i<names.size(); i++) {
                ArrayDouble.D4 varData = new ArrayDouble.D4(1, zDim.getLength(), yDim.getLength(), xDim.getLength());
                double[][][] value = values.get(i);
                boolean[][][] mask = null;
                if (masks != null)
                    mask = masks.get(i);
                for (int xi=0; xi<xDim.getLength(); xi++)
                    for (int yi=0; yi<yDim.getLength(); yi++)
                        for (int zi=0; zi<zDim.getLength(); zi++) {
                            if ((mask == null) || mask[xi][yi][zi])
                                varData.set(0, zi, yi, xi, value[xi][yi][zi]);
                            else
                                varData.set(0, zi, yi, xi, -32768);
                        }
                ncfile.write(names.get(i), dataOrigin, varData);
            }
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    public static void saveProfileVariable(NetcdfFileWriteable ncfile, double time, String name, int dim, double[] values,
                                         boolean writeTime) {
        try {
            ncfile.setRedefineMode(false);
            Dimension tDim = ncfile.findDimension("time"),
                    zDim;
            if (dim==0)
                zDim = ncfile.findDimension("xc");
            else if (dim==1)
                zDim = ncfile.findDimension("yc");
            else zDim = ncfile.findDimension("zc");

            if (writeTime) {
                Array timeData = Array.factory(DataType.DOUBLE, new int[] {1});
                timeData.setDouble(0, time);
                int[] timeOrigin = new int[] {tDim.getLength()};
                ncfile.write("time", timeOrigin, timeData);
            }
            int[] dataOrigin = new int[] {tDim.getLength()-1, 0};

            ArrayDouble.D2 varData = new ArrayDouble.D2(1, zDim.getLength());
            double[] value = values;
            for (int zi=0; zi<zDim.getLength(); zi++) {
                varData.set(0, zi, value[zi]);
            }
            ncfile.write(name, dataOrigin, varData);

        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    public static void saveFishDataVariables(NetcdfFileWriteable ncfile, double time, ArrayList<String> names, ArrayList<double[]> values,
                                         boolean writeTime) {
        try {
            ncfile.setRedefineMode(false);
            Dimension tDim = ncfile.findDimension("time"),
                    xDim = ncfile.findDimension("xc");

            if (writeTime) {
                Array timeData = Array.factory(DataType.DOUBLE, new int[] {1});
                timeData.setDouble(0, time);
                int[] timeOrigin = new int[] {tDim.getLength()};
                ncfile.write("time", timeOrigin, timeData);
            }
            int[] dataOrigin = new int[] {tDim.getLength()-1, 0, 0, 0};

            for (int i=0; i<names.size(); i++) {
                ArrayDouble.D2 varData = new ArrayDouble.D2(1, xDim.getLength());
                double[] value = values.get(i);
                for (int xi=0; xi<xDim.getLength(); xi++)
                    varData.set(0, xi, value[xi]);

                ncfile.write(names.get(i), dataOrigin, varData);
            }
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    public static void saveScalarVariable(NetcdfFileWriteable ncfile, double time, String name, double value, boolean writeTime) {
        ArrayList<String> names = new ArrayList<String>();
        names.add(name);
        ArrayList<Double> values = new ArrayList<Double>();
        values.add(value);
        saveScalarVariables(ncfile, time, names, values, writeTime);
    }

    public static void saveScalarVariables(NetcdfFileWriteable ncfile, double time, ArrayList<String> names, ArrayList<Double> values, boolean writeTime) {
        try {
            ncfile.setRedefineMode(false);
            Dimension tDim = ncfile.findDimension("time");

            if (writeTime) {
                Array timeData = Array.factory(DataType.DOUBLE, new int[] {1});
                timeData.setDouble(0, time);
                int[] timeOrigin = new int[] {tDim.getLength()};
                ncfile.write("time", timeOrigin, timeData);
            }
            int[] dataOrigin = new int[] {tDim.getLength()-1};
            for (int i=0; i<names.size(); i++) {
                ArrayDouble.D1 varData = new ArrayDouble.D1(1);
                double value = values.get(i);
                varData.set(0, value);
                ncfile.write(names.get(i), dataOrigin, varData);
            }
        } catch (Exception ex) {
            ex.printStackTrace();
            System.exit(0);
        }
    }


    public static NetcdfFileWriteable initializeHistogramFile(String f, double[] binEdges, int nCages, String timeUnits) {
        try {
            NetcdfFileWriteable ncfile = NetcdfFileWriteable.createNew(f);
            createTimeDimension(ncfile);
            ncfile.addDimension("bins", binEdges.length);
            ncfile.addDimension("cage", nCages);
            ncfile.create();
            ncfile.setRedefineMode(true);
            createTimeVariable(ncfile, timeUnits);

            // Create variable containing the bin edges:
            ArrayList<Dimension> dimensions = new ArrayList<Dimension>();
            dimensions.add(ncfile.findDimension("bins"));
            Variable beVar = ncfile.addVariable("bin_edges", DataType.DOUBLE, dimensions);
            beVar.addAttribute(new Attribute("scale_factor", 1.0));
            beVar.addAttribute(new Attribute("add_offset", 0.0));

            // Create histogram variable:
            dimensions = new ArrayList<Dimension>();
            dimensions.add(ncfile.findDimension("time"));
            dimensions.add(ncfile.findDimension("cage"));
            dimensions.add(ncfile.findDimension("bins"));
            Variable var = ncfile.addVariable("histogram", DataType.INT, dimensions);
            var.addAttribute(new Attribute("scale_factor", 1.0));
            var.addAttribute(new Attribute("add_offset", 0.0));

            // Write bin edges values:
            ncfile.setRedefineMode(false);
            Array beData = Array.factory(DataType.DOUBLE, new int[] {binEdges.length});
            for (int i=0; i<binEdges.length; i++)
                beData.setDouble(i, binEdges[i]);
            int[] origin = new int[] {0};
            ncfile.write("bin_edges", origin, beData);


            return ncfile;
        } catch (IOException ex) {
            ex.printStackTrace();
            return null;
        } catch (InvalidRangeException e) {
            e.printStackTrace();
            return null;
        }
    }

    public static void saveHistograms(NetcdfFileWriteable ncfile, double time, ArrayList<CageStats> cageStats) {

        try {
            Dimension tDim = ncfile.findDimension("time"),
                    cageDim = ncfile.findDimension("cage"),
                    binsDim = ncfile.findDimension("bins");

            // Write time:
            Array timeData = Array.factory(DataType.DOUBLE, new int[]{1});
            timeData.setDouble(0, time);
            int[] timeOrigin = new int[]{tDim.getLength()};
            ncfile.write("time", timeOrigin, timeData);

            // Write histogram values:
            int[] dataOrigin = new int[] {tDim.getLength()-1, 0, 0};
            ArrayInt.D3 varData = new ArrayInt.D3(1, cageDim.getLength(), binsDim.getLength());
            for (int ci=0; ci<cageDim.getLength(); ci++) {
                int[] histVals = cageStats.get(ci).histogram;
                for (int bi=0; bi< binsDim.getLength(); bi++) {
                    varData.set(0, ci, bi, histVals[bi]);
                }

            }
            ncfile.write("histogram", dataOrigin, varData);

        } catch (InvalidRangeException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
