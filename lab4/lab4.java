import java.util.Arrays;

public class Lab {
    public static void main(String[] args) {
        Program pr = new Program();
        try {
            pr.base();
        } catch (Error e) {
            System.out.println(e.getMessage());
        }
    }
}

class Program {
    double a = -2;
    double b = 2;
    double[] pxk = new double[101];
    double[] f1xk = new double[101];
    double[] f2xk = new double[101];

    Program() {
        for (int i = 0; i < 101; i++) {
            pxk[i] = a + i * (b - a) / 100;
        }
        for (int i = 0; i < 101; i++) {
            f1xk[i] = f1(pxk[i]);
        }
        for (int i = 0; i < 101; i++) {
            f2xk[i] = f2(pxk[i]);
        }
    }

    void base() {
        String s = Arrays.toString(f1xk);
        s = s.substring(1, s.length() - 1);
        System.out.println(s);
        System.out.println("Интерполяционный многочлен в форме Ньютона степени 5 по равностоящим узлам для f1");
        interpolationForF1(5, getUnits(5));
        System.out.println("");
        System.out.println("");
        System.out.println("Интерполяционный многочлен в форме Ньютона степени 10 по равностоящим узлам для f1");
        interpolationForF1(10, getUnits(10));
        System.out.println("");
        System.out.println("");
        System.out.println("Интерполяционный многочлен в форме Ньютона степени 15 по равностоящим узлам для f1");
        interpolationForF1(15, getUnits(15));
        System.out.println("");
        System.out.println("");
        System.out.println("Интерполяционный многочлен в форме Ньютона степени 20 по равностоящим узлам для f1");
        interpolationForF1(20, getUnits(20));
        System.out.println("");
        System.out.println("");
        System.out.println("Интерполяционный многочлен в форме Ньютона степени 5 по чебышевским узлам для f1");
        interpolationForF1(5, getChebUnits(5));
        System.out.println("");
        System.out.println("");
        System.out.println("Интерполяционный многочлен в форме Ньютона степени 10 по чебышевским узлам для f1");
        interpolationForF1(10, getChebUnits(10));
        System.out.println("");
        System.out.println("");
        System.out.println("Интерполяционный многочлен в форме Ньютона степени 15 по чебышевским узлам для f1");
        interpolationForF1(15, getChebUnits(15));
        System.out.println("");
        System.out.println("");
        System.out.println("Интерполяционный многочлен в форме Ньютона степени 20 по чебышевским узлам для f1");
        interpolationForF1(20, getChebUnits(20));
        System.out.println("");
        System.out.println("");
        System.out.println("Интерполяционный многочлен в форме Ньютона степени 5 по равностоящим узлам для f2");
        interpolationForF2(5, getUnits(5));
        System.out.println("");
        System.out.println("");
        System.out.println("Интерполяционный многочлен в форме Ньютона степени 10 по равностоящим узла мдля f2");
        interpolationForF2(10, getUnits(10));
        System.out.println("");
        System.out.println("");
        System.out.println("Интерполяционный многочлен в форме Ньютона степени 15 по равностоящим узлам для f2");
        interpolationForF2(15, getUnits(15));
        System.out.println("");
        System.out.println("");
        System.out.println("Интерполяционный многочлен в форме Ньютона степени 20 по равностоящим узлам для f2");
        interpolationForF2(20, getUnits(20));
        System.out.println("");
        System.out.println("");
        System.out.println("Интерполяционный многочлен в форме Ньютона степени 5 по чебышевским узлам для f2");
        interpolationForF2(5, getChebUnits(5));
        System.out.println("");
        System.out.println("");
        System.out.println("Интерполяционный многочлен в форме Ньютона степени 10 по чебышевским узлам для f2");
        interpolationForF2(10, getChebUnits(10));
        System.out.println("");
        System.out.println("");
        System.out.println("Интерполяционный многочлен в форме Ньютона степени 15 по чебышевским узлам для f2");
        interpolationForF2(15, getChebUnits(15));
        System.out.println("");
        System.out.println("");
        System.out.println("Интерполяционный многочлен в форме Ньютона степени 20 по чебышевским узлам для f2");
        interpolationForF2(20, getChebUnits(20));
        System.out.println("");
    }

    private double[] getUnits(int n) {
        double[] xk = new double[n + 1];
        for (int i = 0; i < n + 1; i++) {
            xk[i] = a + i * (b - a) / n;
        }
        return xk;
    }

    private double[] getChebUnits(int n) {
        double[] cxk = new double[n + 1];
        for (int i = 0; i < n + 1; i++) {
            cxk[i] = (a + b) / 2 + ((b - a) / 2) * Math.cos((2 * i + 1) * Math.PI / (2 * (n + 1)));
        }
        return cxk;
    }

    private double[][] getSeparatedDifferencesForF1(int n, double[] xk) {
        double[][] fy = new double[n + 1][n + 1];
        for (int j = 0; j < n + 1; j++) {
            fy[j][0] = f1(xk[j]);
        }
        int g = 0;
        for (int i = 1; i < n + 1; i++) {
            g++;
            for (int j = 0; j < n + 1 - i; j++) {
                fy[j][i] = (fy[j][i - 1] - fy[j + 1][i - 1]) / (xk[j] - xk[j + g]);
            }
        }
        return fy;
    }

    private double[][] getSeparatedDifferencesForF2(int n, double[] xk) {
        double[][] fy = new double[n + 1][n + 1];
        for (int j = 0; j < n + 1; j++) {
            fy[j][0] = f2(xk[j]);
        }
        int g = 0;
        for (int i = 1; i < n + 1; i++) {
            g++;
            for (int j = 0; j < n + 1 - i; j++) {
                fy[j][i] = (fy[j][i - 1] - fy[j + 1][i - 1]) / (xk[j] - xk[j + g]);
            }
        }
        return fy;
    }

    private double[] interpolationForF1(int n, double[] xk) {
        double pn[] = new double[101];
        double[][] fy = getSeparatedDifferencesForF1(n, xk);
        String s = Double.toString(fy[0][0]);
        String l = "";
        for (int i = 0; i < n; i++) {
            l += "(x - (" + xk[i] + "))";
            s += " + " + l + " * (" + fy[0][i + 1] + ")";
        }
        System.out.println(s);
        for (int i = 0; i < 101; i++) {
            pn[i] = fy[0][0];
            double t = 1;
            for (int j = 0; j < n; j++) {
                t *= pxk[i] - xk[j];
                pn[i] += t * fy[0][j + 1];
            }
        }
        double m = 0;
        for (int i = 0; i < 101; i++) {
            double g = Math.abs(pn[i] - f1xk[i]);
            if (g > m)
                m = g;
        }
        System.out.println("Оценка погрешности " + m);
        return pn;
    }

    private double[] interpolationForF2(int n, double[] xk) {
        double pn[] = new double[101];
        double[][] fy = getSeparatedDifferencesForF2(n, xk);
        String s = Double.toString(fy[0][0]);
        StringBuilder l = new StringBuilder();
        for (int i = 0; i < n; i++) {
            if(i < 1)
                l.append("(x - (").append(xk[i]).append("))");
            else
                l.append("*(x - (").append(xk[i]).append("))");
            s += " + " + l + " * (" + fy[0][i + 1] + ")";
        }
        System.out.println(s);
        for (int i = 0; i < 101; i++) {
            pn[i] = fy[0][0];
            double t = 1;
            for (int j = 0; j < n; j++) {
                t *= pxk[i] - xk[j];
                pn[i] += t * fy[0][j + 1];
            }
        }
        double m = 0;
        for (int i = 0; i < 101; i++) {
            double g = Math.abs(pn[i] - f2xk[i]);
            if (g > m)
                m = g;
        }
        System.out.println("Оценка погрешности " + m);
        return pn;
    }

    private double f1(double x) {
        return Math.sin(Math.cos(x));
    }

    private double f2(double x) {
        return Math.abs(Math.abs(x) - 1);
    }
}
