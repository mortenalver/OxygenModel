package fishmodel;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Reads all simulation configuration from a JSON file.
 * No external JSON library is required: a minimal parser is embedded below.
 *
 * Usage:
 *   SimConfig cfg = new SimConfig("path/to/loenngrunnen.json");
 *   // then use cfg.dxy, cfg.inputDataDir, cfg.cagePositions, …
 */
public class SimConfig {

    // ---- File paths -------------------------------------------------------
    /** Base directory for field-measurement input files (nc, csv). */
    public final String inputDataDir;
    /** NetCDF file name for current/O2 field observations. */
    public final String inDataFileName;
    /** NetCDF file name for inventory (fish weight, count, feed) data. */
    public final String invDataFileName;
    /** Directory containing the spatially variable current-magic fields. */
    public final String currentMagicDir;
    /** File name of the current-magic field NetCDF. */
    public final String currentMagicFileName;
    /** Output directory for simulation result files. */
    public final String saveDir;
    /** Optional postfix appended to every output file name. */
    public final String simNamePostfix;
    /** Short scenario label used in file/folder names. */
    public final String outputPrefix;

    // ---- Mode flags -------------------------------------------------------
    public final boolean useFieldInputs;
    public final boolean useCurrentMagic;
    public final boolean useNewConsumptionModel;
    public final boolean useDigestiveModel;
    /** If true, ignore inventory and use fixed biomass/feeding (Fluent validation). */
    public final boolean testEastCurrAgainstFluent;
    public final boolean cage6BiomassOnly;
    public final boolean dynamicCurrReductionFactor;
    public final boolean includeHypoxiaAvoidance;
    /** Reduce ambient O2 and feed near the edges of the farm domain. */
    public final boolean varyAmbient;
    /** Use depth-dependent affinity profile (true) or uniform (false). */
    public final boolean useVerticalDist;
    public final boolean maskO2WhenSaving;
    public final boolean storeO2Histograms;

    // ---- Numeric parameters ----------------------------------------------
    /** Horizontal model resolution (m). */
    public final double dxy;
    /** Vertical model resolution (m). Defaults to dxy when omitted. */
    public final double dz;
    /** Number of daily simulation periods. */
    public final int nSim;
    /** Skip this many simulations at the start (but count them in numbering). */
    public final int startAt;
    public final double extraConsumptionMultiplier;
    public final double currentReductionFactor;
    /** Fraction of fish that do NOT follow the feed distribution. */
    public final double o2EvenFraction;
    /** Initial external current speed (m/s). */
    public final double currentSpeedInit;
    /** Water temperature (°C). */
    public final double T_w;
    /** Ambient dissolved oxygen (mg/L). */
    public final double avO2;

    // ---- Simulation start date -------------------------------------------
    public final int startYear;
    /** Start month, 1-based (January = 1). */
    public final int startMonth;
    public final int startDay;

    // ---- Farm geometry ---------------------------------------------------
    public final int cageGridRows;
    public final int cageGridCols;
    /** [nCages][2] — each row is {gridRow, gridCol} for one cage. */
    public final int[][] cagePositions;
    public final double farmRotation;
    public final double frameSize;
    public final double outerPadding;
    /** [surfaceRadius, bottomRadius] — spaghetti cage radii (m). */
    public final double[] spaghettiRad;
    /** [bottomTauDepth, bottomNetDepth] — spaghetti cage depths (m). */
    public final double[] spaghettiDepth;
    public final double depthDomain;
    public final double fishMaxDepth;
    /** Feeding start time — seconds from midnight. */
    public final int feedStartSeconds;
    /** Feeding end time — seconds from midnight. */
    public final int feedEndSeconds;

    // ---- Sensor layout ---------------------------------------------------
    public final double[] sensorDepths;
    /** Angles (degrees from north) of horizontal sensor positions. */
    public final double[] sensorAngles;
    /** Distance of sensors from cage centre (m). */
    public final double sensorRad;

    // ---- Fish affinity vertical profile ----------------------------------
    /** Depths (m) at which the affinity profile values are given. */
    public final double[] affinityDepths;
    /** Relative affinity values at {@link #affinityDepths}. */
    public final double[] affinityValues;

    // ---- Output settings -------------------------------------------------
    /** How often to store feed-distribution data (s). */
    public final int storeIntervalFeed;
    /** How often to store general info (s). */
    public final int storeIntervalInfo;
    /** How often to re-evaluate O2-avoidance affinity (time steps). */
    public final int checkAvoidanceInterval;

    // ======================================================================
    // Constructor
    // ======================================================================

    public SimConfig(String configPath) throws IOException {
        String json = new String(Files.readAllBytes(Paths.get(configPath)), StandardCharsets.UTF_8);
        Map<String, Object> m = parseObject(new Tok(json));

        inputDataDir              = str(m, "inputDataDir",         "");
        inDataFileName            = str(m, "inDataFileName",        "loenngrunnen_data_mid_fix.nc");
        invDataFileName           = str(m, "invDataFileName",       "Loenngrunnen_inv.nc");
        currentMagicDir           = str(m, "currentMagicDir",       "");
        currentMagicFileName      = str(m, "currentMagicFileName",  "lti_currents_21may-24jun.nc");
        saveDir                   = str(m, "saveDir",               "./output/");
        simNamePostfix            = str(m, "simNamePostfix",        "");
        outputPrefix = str(m, "outputPrefix",           "sim");

        useFieldInputs            = bool(m, "useFieldInputs",            true);
        useCurrentMagic           = bool(m, "useCurrentMagic",           false);
        useNewConsumptionModel    = bool(m, "useNewConsumptionModel",    true);
        useDigestiveModel         = bool(m, "useDigestiveModel",         false);
        testEastCurrAgainstFluent = bool(m, "testEastCurrAgainstFluent", false);
        cage6BiomassOnly          = bool(m, "cage6BiomassOnly",          false);
        dynamicCurrReductionFactor= bool(m, "dynamicCurrReductionFactor",false);
        includeHypoxiaAvoidance   = bool(m, "includeHypoxiaAvoidance",   true);
        varyAmbient               = bool(m, "varyAmbient",               false);
        useVerticalDist           = bool(m, "useVerticalDist",           true);
        maskO2WhenSaving          = bool(m, "maskO2WhenSaving",          false);
        storeO2Histograms         = bool(m, "storeO2Histograms",         true);

        dxy                       = num(m, "dxy",                     4.0);
        dz                        = num(m, "dz",                      dxy);
        nSim                      = (int) num(m, "nSim",              49);
        startAt                   = (int) num(m, "startAt",           0);
        extraConsumptionMultiplier= num(m, "extraConsumptionMultiplier", 1.0);
        currentReductionFactor    = num(m, "currentReductionFactor",  0.8);
        o2EvenFraction            = num(m, "o2EvenFraction",          0.75);
        currentSpeedInit          = num(m, "currentSpeedInit",        0.08);
        T_w                       = num(m, "T_w",                     14.0);
        avO2                      = num(m, "avO2",                    9.0);

        startYear                 = (int) num(m, "startYear",         2025);
        startMonth                = (int) num(m, "startMonth",        5);
        startDay                  = (int) num(m, "startDay",          21);

        cageGridRows              = (int) num(m, "cageGridRows",      6);
        cageGridCols              = (int) num(m, "cageGridCols",      2);
        farmRotation              = num(m, "farmRotation",            -84.35);
        frameSize                 = num(m, "frameSize",               110.0);
        outerPadding              = num(m, "outerPadding",            75.0);
        depthDomain               = num(m, "depthDomain",             38.0);
        fishMaxDepth              = num(m, "fishMaxDepth",            38.0);
        feedStartSeconds          = (int) num(m, "feedStartSeconds",  7 * 3600);
        feedEndSeconds            = (int) num(m, "feedEndSeconds",    19 * 3600);
        sensorRad                 = num(m, "sensorRad",               10.0);

        storeIntervalFeed         = (int) num(m, "storeIntervalFeed",       7200);
        storeIntervalInfo         = (int) num(m, "storeIntervalInfo",         60);
        checkAvoidanceInterval    = (int) num(m, "checkAvoidanceInterval",    30);

        spaghettiRad   = doubleArr(m, "spaghettiRad",    new double[]{30.0, 30.0});
        spaghettiDepth = doubleArr(m, "spaghettiDepth",  new double[]{20.0, 38.0});
        sensorDepths   = doubleArr(m, "sensorDepths",    new double[]{5.0, 10.0, 15.0});
        sensorAngles   = doubleArr(m, "sensorAngles",    new double[]{45.0});
        affinityDepths = doubleArr(m, "affinityDepths",  new double[]{0.5, 27.5});
        affinityValues = doubleArr(m, "affinityValues",  new double[]{1.0, 1.0});
        cagePositions  = intPairArr(m, "cagePositions",  new int[0][]);
    }

    // ======================================================================
    // Typed accessors for the parsed map
    // ======================================================================

    private static String str(Map<String,Object> m, String k, String def) {
        Object v = m.get(k); return v instanceof String ? (String)v : def;
    }
    private static boolean bool(Map<String,Object> m, String k, boolean def) {
        Object v = m.get(k); return v instanceof Boolean ? (Boolean)v : def;
    }
    private static double num(Map<String,Object> m, String k, double def) {
        Object v = m.get(k); return v instanceof Double ? (Double)v : def;
    }

    @SuppressWarnings("unchecked")
    private static double[] doubleArr(Map<String,Object> m, String k, double[] def) {
        Object v = m.get(k);
        if (!(v instanceof List)) return def;
        List<Object> list = (List<Object>)v;
        double[] a = new double[list.size()];
        for (int i = 0; i < a.length; i++) a[i] = (Double)list.get(i);
        return a;
    }

    @SuppressWarnings("unchecked")
    private static int[][] intPairArr(Map<String,Object> m, String k, int[][] def) {
        Object v = m.get(k);
        if (!(v instanceof List)) return def;
        List<Object> outer = (List<Object>)v;
        int[][] a = new int[outer.size()][2];
        for (int i = 0; i < outer.size(); i++) {
            List<Object> row = (List<Object>)outer.get(i);
            a[i][0] = ((Double)row.get(0)).intValue();
            a[i][1] = ((Double)row.get(1)).intValue();
        }
        return a;
    }

    // ======================================================================
    // Minimal recursive-descent JSON parser (no external dependencies)
    // ======================================================================

    private static Map<String,Object> parseObject(Tok t) {
        t.expect('{');
        Map<String,Object> m = new LinkedHashMap<>();
        t.skipWs();
        if (t.peek() == '}') { t.advance(); return m; }
        while (true) {
            t.skipWs();
            String key = t.readString();
            t.skipWs(); t.expect(':'); t.skipWs();
            m.put(key, parseValue(t));
            t.skipWs();
            char c = t.advance();
            if (c == '}') break;
            if (c != ',') throw new RuntimeException("JSON: expected ',' or '}', got '" + c + "' near pos " + t.pos);
        }
        return m;
    }

    private static Object parseValue(Tok t) {
        t.skipWs();
        char c = t.peek();
        if (c == '"')  return t.readString();
        if (c == '{')  return parseObject(t);
        if (c == '[')  return parseArray(t);
        if (c == 't')  { t.expect("true");  return Boolean.TRUE; }
        if (c == 'f')  { t.expect("false"); return Boolean.FALSE; }
        if (c == 'n')  { t.expect("null");  return null; }
        return t.readNumber();
    }

    private static List<Object> parseArray(Tok t) {
        t.expect('[');
        List<Object> list = new ArrayList<>();
        t.skipWs();
        if (t.peek() == ']') { t.advance(); return list; }
        while (true) {
            t.skipWs();
            list.add(parseValue(t));
            t.skipWs();
            char c = t.advance();
            if (c == ']') break;
            if (c != ',') throw new RuntimeException("JSON: expected ',' or ']', got '" + c + "' near pos " + t.pos);
        }
        return list;
    }

    /** Minimal tokeniser over the JSON string. */
    private static class Tok {
        final String src;
        int pos;
        Tok(String src) { this.src = src; }

        void skipWs() {
            while (pos < src.length() && Character.isWhitespace(src.charAt(pos))) pos++;
        }
        char peek() { skipWs(); return pos < src.length() ? src.charAt(pos) : 0; }
        char advance() { skipWs(); return pos < src.length() ? src.charAt(pos++) : 0; }

        void expect(char c) {
            skipWs();
            if (pos >= src.length() || src.charAt(pos) != c)
                throw new RuntimeException("JSON: expected '" + c + "' at pos " + pos);
            pos++;
        }
        void expect(String s) {
            skipWs();
            for (int i = 0; i < s.length(); i++) {
                if (pos >= src.length() || src.charAt(pos) != s.charAt(i))
                    throw new RuntimeException("JSON: expected literal '" + s + "' at pos " + (pos - i));
                pos++;
            }
        }

        String readString() {
            expect('"');
            StringBuilder sb = new StringBuilder();
            while (pos < src.length()) {
                char c = src.charAt(pos++);
                if (c == '"') return sb.toString();
                if (c == '\\' && pos < src.length()) {
                    char e = src.charAt(pos++);
                    switch (e) {
                        case '"':  sb.append('"');  break;
                        case '\\': sb.append('\\'); break;
                        case '/':  sb.append('/');  break;
                        case 'n':  sb.append('\n'); break;
                        case 'r':  sb.append('\r'); break;
                        case 't':  sb.append('\t'); break;
                        default:   sb.append(e);
                    }
                } else {
                    sb.append(c);
                }
            }
            throw new RuntimeException("JSON: unterminated string");
        }

        Double readNumber() {
            skipWs();
            int start = pos;
            if (pos < src.length() && src.charAt(pos) == '-') pos++;
            while (pos < src.length()) {
                char c = src.charAt(pos);
                if (Character.isDigit(c) || c == '.' || c == 'e' || c == 'E' || c == '+' || c == '-')
                    pos++;
                else break;
            }
            if (pos == start) throw new RuntimeException("JSON: expected number at pos " + pos);
            return Double.parseDouble(src.substring(start, pos));
        }
    }
}
