package fishmodel.pellets;

import save.SaveNetCDF;
import ucar.nc2.NetcdfFileWriteable;

import java.io.IOException;

/**
 * Model for horizontal distribution of feed using rotating spreader. (CREATE project with NOFIMA)
 */
public class PelletSpreaderModel {

    static double[][] param = new double[][] {
            {1.2904, 0.0059, 0.0424, 1.0628, 0.0185, 0.0639, 0.9015, 1.0189, 0.2577},
            {0.7254, 0.1307, 0.0972, 1.6763, 0.1304, 0.1476, 0.7822, 0.5951, 0.7905}
    };

    static final double AVERAGE_PELLET_SPEED = 12; // m/s average speed from pellets leave the spreader until they hit the water

    static final double FORCE_PROPORTIONALITY = 0.5*1;

    static final double SPREAD_FACTOR = 0.25; // The amount of spread of wind-blown feed per cell length the feed is displaced

    public static void setPelletDist(double[][] feedingDist, int centerX, int centerY, double dxy, double feederAngle,
                                     double airSpeed, int spreaderType, boolean tiltUp, double[] windSpeed, boolean[][][] mask) {

        double sum = 0;

        for (int i=0; i<feedingDist.length; i++)
            for (int j=0; j<feedingDist[0].length; j++) {
                feedingDist[i][j] = 0;
            }

        final int PAD = 10;
        double[][] tmp = new double[feedingDist.length+2*PAD][feedingDist[0].length+2*PAD];


        for (int i=0; i<feedingDist.length; i++)
            for (int j=0; j<feedingDist[i].length; j++) {
                
                double distX = dxy*(i-centerX), distY = dxy*(j-centerY);
                double distance = Math.sqrt(distX*distX + distY*distY);
                double angle = (180./Math.PI)*Math.atan2(distX, distY) - feederAngle;
                while (angle < 0) angle += 360;

                // Effect of wind. We assume the wind causes a constant force on each pellet in the wind direction (relative
                // to the case of no wind, for which the spreader model is valid).
                // This gives a delta distance travelled equal to 0.5*(F/m)*t^2.
                // Assume that the time for each pellet to land is proportional to its original distance to travel.
                // The acceleration F/m is assumed to be proportional to wind speed.
                // We can calculate the delta separately along the two wind components:
                double airTime = distance/AVERAGE_PELLET_SPEED;
                double airTimeSquared = airTime*airTime;
                // Calculate wind driven delta movement in cells:
                double[] windDelta = new double[] {windSpeed[0]*FORCE_PROPORTIONALITY*airTimeSquared/dxy,
                        windSpeed[1]*FORCE_PROPORTIONALITY*airTimeSquared/dxy};
                //System.out.println("airTime = "+airTime);
                //System.out.println("windDelta = "+windDelta[0]+" , "+windDelta[1]);

                // Amount of feed to add in this iteration:
                double amount = distance > 0 ? (1.0/distance)*pelletSpread(distance, angle, airSpeed, spreaderType, tiltUp) : 0;

                double windDeltaAbs = Math.sqrt(windDelta[0]*windDelta[0] + windDelta[1]*windDelta[1]);

                if (windDeltaAbs == 0) { // No wind, simple addition
                    feedingDist[i][j] += amount;
                }
                else {
                    double new_i = i+windDelta[0];
                    double new_j = j+windDelta[1];
                    double spread = 1 + SPREAD_FACTOR*windDeltaAbs;

                    // Break up the wind delta in a distribution between two cell numbers:
                    /*int i_low = (int)Math.floor(windDelta[0]);
                    double i_fraction = windDelta[0] - Math.floor(windDelta[0]);
                    int j_low = (int)Math.floor(windDelta[1]);
                    double j_fraction = windDelta[1] - Math.floor(windDelta[1]);
                    */
                    double hsum = 0;
                    for (int ii=0; ii<tmp.length; ii++)
                        for (int jj=0; jj<tmp[0].length; jj++) {
                            double dH = Math.sqrt(((double)(ii-PAD)-new_i)*((double)(ii-PAD)-new_i) + ((double)(jj-PAD)-new_j)*((double)(jj-PAD)-new_j));
                            tmp[ii][jj] = Math.max(0, 1-dH/spread);
                            hsum += tmp[ii][jj];
                        }

                    if (hsum > 0)
                        for (int ii=0; ii<feedingDist.length; ii++)
                            for (int jj=0; jj<feedingDist[0].length; jj++) {
                                feedingDist[ii][jj] += tmp[ii+PAD][jj+PAD]*amount/hsum;
                            }

                            //i_low = 0;
                    //j_low = 0;

                    //if ((i!= 14) || (j != 14)) amount = 0;

                    // Spread among four cells:
                    /*if (i_fraction < 0.1) i_fraction = 0;
                    if (i_fraction > 0.9) i_fraction = 1;
                    if (j_fraction < 0.1) j_fraction = 0;
                    if (j_fraction > 0.9) j_fraction = 1;*/

                    /*if (j == 60)
                        System.out.println("i = "+i+"; actual i = "+(i+windDelta[0])+"; i_fraction = "+i_fraction);

                    //addIfInside(feedingDist, mask, i+i_low, j+j_low, amount);
                    addIfInside(feedingDist, mask, i+i_low, j+j_low, amount*(1-i_fraction)*(1-j_fraction));
                    addIfInside(feedingDist, mask, i+i_low+1, j+j_low, amount*i_fraction*(1-j_fraction));
                    addIfInside(feedingDist, mask, i+i_low, j+j_low+1, amount*(1-i_fraction)*j_fraction);
                    addIfInside(feedingDist, mask, i+i_low+1, j+j_low+1, amount*i_fraction*j_fraction);
                     */

                    /*if ((mask == null) || mask[i][j][0]) { // With masking, only feed inside the cage/tank:
                        feedingDist[i][j] = distance > 0 ? amount : 0;
                        sum += feedingDist[i][j];
                    }*/
                }
            }

            for (int i=0; i<feedingDist.length; i++)
                for (int j=0; j<feedingDist[0].length; j++) {
                    sum += feedingDist[i][j];
                }
            //double sum2 = 0;
            for (int i=0; i<feedingDist.length; i++)
                for (int j=0; j<feedingDist[0].length; j++) {
                    feedingDist[i][j] /= sum;
                    //sum2 += feedingDist[i][j];
                }

    }

    private static void addIfInside(double[][] array, boolean[][][] mask, int i, int j, double value) {
        if ((i >= 0) && (j >= 0) && (i < array.length) && (j < array[0].length))
            if ((mask == null) || mask[i][j][0]) { // With masking, only feed inside the cage/tank:
                array[i][j] += value;
            }
    }

    private static double pelletSpread(double distance, double angle, double airspeed, int spreaderType, boolean tiltUp) {
        distance = Math.pow(distance, param[spreaderType][8]);

        double halfangle;
        if (angle < 180)
            halfangle = angle/180;
        else
            halfangle = (360 - angle)/180;

        double centerFwd = param[spreaderType][0] + airspeed*param[spreaderType][1];
        double centerBack = param[spreaderType][3] + airspeed*param[spreaderType][4];
        if (tiltUp) {
            centerFwd = centerFwd*param[spreaderType][6];
            centerBack = centerBack*param[spreaderType][7];
        }
        //centerBack = centerFwd;
        double spreadFwd = centerFwd*param[spreaderType][2];
        double spreadBack = centerBack*param[spreaderType][5];
        //spreadBack = spreadFwd;
        double center = centerBack*halfangle + centerFwd*(1-halfangle);
        double spread = spreadBack*halfangle + spreadFwd*(1-halfangle);

        return (1/spread*Math.sqrt(2*Math.PI))*Math.exp(-0.5 * ((distance - center) / spread) * ((distance - center) / spread));
    }

    public static void copy(double[][] from, double[][][] to) {
        for (int i=0; i<from.length; i++)
            for (int j=0; j<from[i].length; j++)
                to[i][j][0] = from[i][j];
    }

    public static void main(String[] args) {
        int[] cageDims = new int[] {120, 120, 1};
        double dxy = 0.2/1;
        double[] windSpeed = new double[] {0, 0};
        NetcdfFileWriteable ncfile = SaveNetCDF.initializeFile("C:\\temp\\testfile_tiltdown2.nc", cageDims, 1, 1);
        SaveNetCDF.createCageVariables(ncfile, "feed");
        double[][] feedingDist = new double[cageDims[0]][cageDims[1]];
        double[][][] feed = new double[cageDims[0]][cageDims[1]][cageDims[2]];

        double ws = 0, wd = 0;


        for (int i=0; i<1/*200*/; i++) {
            System.out.println("i = "+i);

            // Update wind speed and direction at random:
            ws = 0*ws + 13*(2.5+Math.sin(((double)i) / 5));//(Math.random()-0.5);
            wd = wd + 0.0285;//0.25*(Math.random()-0.5);
            if (i < 10)
                ws = ws*((double)i/10.);

            windSpeed[0] = ws*Math.cos(wd);
            windSpeed[1] = ws*Math.sin(wd);
            System.out.println(windSpeed[0]+" : "+windSpeed[1]);
            setPelletDist(feedingDist, cageDims[0]/2, cageDims[1]/2, dxy, 0, 25, 0, false, windSpeed, null);
            copy(feedingDist, feed);
            SaveNetCDF.saveCageVariable(ncfile, i, "feed", feed, null, true);
        }

        try {
            ncfile.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
