import static java.lang.Math.pow;

public class Lab {
    public static void main(String[] args) {
        Program pr = new Program();
        try {
            pr.base();
        }catch (Error e){
            System.out.println(e.getMessage());
        }
    }
}
class Program {
    void base() {
        double e  = 1.0E-7;
        double x = methodOfDichotomy(0.1);
        System.out.println();
        System.out.println();
        methodOfSimpleIteration(e, x, true);
        System.out.println();
        System.out.println();
        methodOfNewton(e, x);
        System.out.println();
        System.out.println();
        methodOfStevenson(e, x);
    }

    private  double f(double x) {
        return Math.exp(x)-pow(x,3)+3*pow(x,2) - 2 * x - 3;
    }

    private double df(double x) {
        return Math.exp(x)-3 * pow(x,2)+6*x - 2;
    }

    private  double sf(double x) {
        return -Math.sqrt((-Math.exp(x) + pow(x,3) + 2 * x + 3)/3.0);
    }

    private  double kvsf(double x) {
        return (-Math.exp(x) + pow(x,3) + 2 * x + 3)/3.0;
    }

    private double methodOfDichotomy(double e) {
        double a = -1;
        double b = 0.01;
        double x = 0;
        int k = 0;
        System.out.println("Метод дихотомии: ");
        System.out.print("Номер итерации: " + k + "        ");
        System.out.printf("a(k) = %.1f" + "       ", a);
        System.out.printf("b(k) = %.2f" + "       ", b);
        System.out.printf("f(a(k)) = %.8f" + "       ", f(a));
        System.out.printf("f(b(k)) = %.8f" + "       ", f(b));
        System.out.printf("((a(k) + b(k)) / 2 = %.4f" + "       ", (a+b)/2);
        System.out.printf("b(k) - a(k) = %.2f\n",  b - a);
        while (Math.abs(a-b) > e) {
            x = (a+b)/2;
            if (f(a) * f(x)<= 0)
                b = x;
            else
                a = x;
            k++;
            System.out.print("Номер итерации: " + k + "        ");
            System.out.printf("a(k) = %.6f" + "       ", a);
            System.out.printf("b(k) = %.6f" + "       ", b);
            System.out.printf("f(a(k)) = %.8f" + "       ", f(a));
            System.out.printf("f(b(k)) = %.8f" + "       ", f(b));
            System.out.printf("((a(k) + b(k)) / 2 = %.6f" + "       ", x);
            System.out.printf("b(k) - a(k) = %.8f\n",  b - a);
        }
        return x;
    }

    private void methodOfSimpleIteration(double e, double x, boolean b){
        int k = 0;
        double xk = Double.MAX_VALUE;
        System.out.println("Метод простой итерации");
        System.out.printf("Номер итерации: " + k + "        ");
        System.out.printf("x(k) = %.8f\n", x);
        while ( Math.abs(xk-x) > e){
            k++;
            xk = x;
            x = sf(x);
            if(b){
                System.out.print("Номер итерации: " + k + "        ");
                System.out.printf("x(k) = %.16f" + "       ", x);
                System.out.printf("|x(k) - x(k-1)| = %.16f\n",  Math.abs(xk-x));
            }
            if(k>20000)
                throw new Error("Превышен лимит итераций, метод скорее всего не сходится");
        }
    }

    private void methodOfNewton(double e, double x){
        int k = 0;
        double xk = Double.MAX_VALUE;
        System.out.println("Метод Ньютона");
        System.out.printf("Номер итерации: " + k + "        ");
        System.out.printf("x(k) = %.8f\n", x);
        while ( Math.abs(xk-x) > e){
            k++;
           if (df(x) == 0) {
                break;
            }
            xk = x;
            x -= f(x) / df(x);
            System.out.print("Номер итерации: " + k + "        ");
            System.out.printf("x(k) = %.16f" + "       ", x);
            System.out.printf("|x(k) - x(k-1)| = %.16f\n",  Math.abs(xk-x));
        }
    }

    private void methodOfStevenson(double e, double x) {
        int k = 0;
        double xk = Double.MAX_VALUE;
        System.out.println("Метод Стеффенсена");
        System.out.printf("Номер итерации: " + k + "        ");
        System.out.printf("x(k) = %.8f\n", x);
        while (Math.abs(xk - x) > e) {
            k++;
            xk = x;
            x = (x * sf(sf(x)) - kvsf(x)) / (sf(sf(x)) - 2 * sf(x) + x);
            System.out.print("Номер итерации: " + k + "        ");
            System.out.printf("x(k) = %.16f" + "       ", x);
            System.out.printf("|x(k) - x(k-1)| = %.16f\n", Math.abs(xk - x));
        }
    }
}
