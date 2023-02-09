package fishmodel.pellets;

/**
 * Implementation of formula from Garcia & Gordon (1992)
 */
public class OxygenSolubility {


    /**
     * Calculate the solubility of oxygen in 35ppt sea water as a function of temperature
     * @param T The temperature in C
     * @return The solubility in mg/l
     */
    public static double getOxygenSolubility(double T) {
        double T_s = Math.log((298.15 - T)/(273.15 + T));
        double S = 35;
        double lnC = 5.80818 + 3.20684*T_s + 4.11890*T_s*T_s +
                4.93845*(0*T_s*T_s + T_s*T_s*T_s) + 1.0567*Math.pow(T_s, 4) +
                1.41575*Math.pow(T_s, 5)
                + S*(-7.01211e-3 -7.25958e-3*T_s - 7.93334e-3*T_s*T_s
                    - 5.54491e-3*T_s*T_s*T_s) - 1.32412e-7*S*S;
        return 0.001*32.0*1.035*Math.exp(lnC);
    }

    public static void main(String[] args) {
        for (int i=0; i<20; i++) {
            System.out.println("Temp "+(i*1.0)+": sol = "+getOxygenSolubility(i*1.0));
        }
    }
}
