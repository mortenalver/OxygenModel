package fishmodel.enkf;

import mpi.*;

public class MpiHandler {

    int rank, N;

    public MpiHandler(String[] args) {
        MPI.Init(args);
        rank = MPI.COMM_WORLD.Rank();
        N = MPI.COMM_WORLD.Size();
    }

    public int getRank() {
        return rank;
    }

    public int getN() {
        return N;
    }

    public double[][] gatherStateToRank0(double[][][] field, double[] parVal) {
        int nPar = parVal.length;
        if (rank==0) {
            double[][] allStates = new double[N][];
            //double[][] allParVals = null;
            //if (nPar > 0)
            //    allParVals = new double[N][nPar];
            double[] vector = Util.reshapeArray(field);
            allStates[0] = new double[vector.length + nPar];
            System.arraycopy(vector, 0, allStates[0], 0, vector.length);
            if (nPar > 0)
                System.arraycopy(parVal, 0, allStates[0], vector.length, nPar);
            for (int i=1; i<N; i++) {
                MPI.COMM_WORLD.Recv(vector, 0, vector.length, MPI.DOUBLE, i, 1);
                allStates[i] = new double[vector.length + nPar];
                System.arraycopy(vector, 0, allStates[i], 0, vector.length);
                if (nPar > 0) {
                    double[] parValRecv = new double[nPar];
                    MPI.COMM_WORLD.Recv(parValRecv, 0, nPar, MPI.DOUBLE, i, 2);
                    System.arraycopy(parValRecv, 0, allStates[i], vector.length, nPar);
                }
                //System.out.println("Received data from rank "+i);
            }

            /*for (int i=0; i< allStates.length; i++) {
                for (int j=0; j<allStates[i].length; j++) {
                    System.out.print(allStates[i][j]+"\t");
                }
                System.out.println("");
            }*/


            return allStates;
        }
        else {
            double[] vector = Util.reshapeArray(field);
            MPI.COMM_WORLD.Send(vector, 0, vector.length, MPI.DOUBLE, 0, 1);
            if (nPar > 0)
                MPI.COMM_WORLD.Send(parVal, 0, parVal.length, MPI.DOUBLE, 0, 2);
            return null;
        }
    }

    public void distributeAnalysisFromRank0(double[][] X_a, double[][][] field, double[] parVal, int[] dims, int nPar) {
        //System.out.println("DistributeAnalysis: "+X_a.length+" x "+X_a[0].length);
        // Send out all analysis vectors except the one for rank 0:
        for (int i=1; i<X_a.length; i++) {
            MPI.COMM_WORLD.Send(X_a[i], 0, X_a[i].length, MPI.DOUBLE, i, 1);
        }
        Util.reshapeInto3d(X_a[0], field, dims);
        if (nPar > 0) {
            int numel = dims[0]*dims[1]*dims[2];
            for (int i=0; i<nPar; i++)
                parVal[i] = X_a[0][numel + i];
        }
    }

    public void receiveAnalysisFromRank0(double[][][] field, double[] parVal, int[] dims, int nPar) {
        int numel = dims[0]*dims[1]*dims[2];
        int nValues = numel + nPar;
        double[] X_a = new double[nValues];
        MPI.COMM_WORLD.Recv(X_a, 0, X_a.length, MPI.DOUBLE, 0, 1);
        Util.reshapeInto3d(X_a, field, dims);
        if (nPar > 0) {
            for (int i=0; i<nPar; i++)
                parVal[i] = X_a[numel + i];
        }
    }
}
