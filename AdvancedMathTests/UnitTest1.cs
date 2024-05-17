using AdvancedMath;
using System.Globalization;

namespace MathLibrary.Tests
{
    [TestFixture]
    public class MuavraLaplaceSolverWithStableKTests
    {
        [Test]
        public void CalculateProbability_Returns_CorrectProbability()
        {
            // Arrange
            double p = 0.8;
            int n = 100;
            int k = 75;
            var solver = new MuavraLaplaceSolverWithStableK(p, k, n);

            // Act
            double probability = solver.CalculateProbability();

            // Assert
            Assert.AreEqual(0.043746, probability, 0.005);
        }
    }

    [TestFixture]
    public class ArrayOperationsTests
    {
        [Test]
        public void CalculateArrayMean_Returns_CorrectMean()
        {
            // Arrange
            double[] values = { 2, 4, 6, 8, 10 };

            // Act
            double mean = ArrayOperations.CalculateArrayMean(values);

            // Assert
            Assert.AreEqual(6, mean);
        }

        [Test]
        public void CalculateStandardArrayDeviation_Returns_CorrectDeviation()
        {
            // Arrange
            double[] values = { 2, 4, 6, 8, 10 };

            // Act
            double deviation = ArrayOperations.CalculateStandardArrayDeviation(values);

            // Assert
            Assert.AreEqual(2.828427, deviation, 0.000001);
        }
    }

    [TestFixture]
    public class ExtraMathOperationsTests
    {
        [Test]
        public void NewtonRaphsonMethod_Returns_CorrectRoot()
        {
            // Arrange
            Func<double, double> function = x => Math.Pow(x, 3) - 2 * x - 5;
            Func<double, double> derivative = x => 3 * Math.Pow(x, 2) - 2;
            double initialGuess = 2;
            double tolerance = 0.0001;
            int maxIterations = 100;

            // Act
            double root = ExtraMathOperations.NewtonRaphsonMethod(function, derivative, initialGuess, tolerance, maxIterations);

            // Assert
            Assert.AreEqual(2.09455, root, 0.00001);
        }

        [Test]
        public void NumericalIntegration_Returns_CorrectIntegral()
        {
            // Arrange
            Func<double, double> squareFunction = x => x * x;
            double lowerBound = 0;
            double upperBound = 1;
            int numIntervals = 100;

            // Act
            double integral = ExtraMathOperations.NumericalIntegration(squareFunction, lowerBound, upperBound, numIntervals);

            // Assert
            Assert.AreEqual(0.333333, integral, 0.0001);
        }
    }

    [TestFixture]
    public class DifferentiationTests
    {
        [Test]
        public void RungeKuttaSolver_Returns_CorrectSolution()
        {
            // Arrange
            Func<double, double, double> differentialEquation = (x, y) => x * y;
            double initialValueX = 0;
            double initialValueY = 1;
            double stepSize = 0.1;
            double endTime = 1;

            // Act
            double[] resultY = Differentiaton.RungeKuttaSolver(differentialEquation, initialValueX, initialValueY, stepSize, endTime);

            // Assert
            Assert.AreEqual(1.645, resultY[resultY.Length - 1], 0.005);
            Assert.AreEqual(1.375, resultY[resultY.Length - 3], 0.005);
            Assert.AreEqual(1.195, resultY[resultY.Length - 5], 0.005);
        }

        [Test]
        public void SolveDifferentialEquationEuler_Returns_CorrectSolution()
        {
            // Arrange
            Func<double, double, double> derivative = (x, y) => x * y;
            double initialValueX = 0;
            double initialValueY = 1;
            double stepSize = 0.1;
            double endTime = 1;

            // Act
            double[] resultY = Differentiaton.SolveDifferentialEquationEuler(derivative, initialValueX, initialValueY, stepSize, endTime);

            // Assert
            Assert.AreEqual(1.55, resultY[resultY.Length - 1], 0.01);
            Assert.AreEqual(1.31, resultY[resultY.Length - 3], 0.01);
            Assert.AreEqual(1.15, resultY[resultY.Length - 5], 0.01);
        }
    }

    [TestFixture]
    public class ComplexNumberTests
    {
        [Test]
        public void ComplexNumber_Addition_Returns_CorrectResult()
        {
            // Arrange
            ComplexNumber c1 = new ComplexNumber(1, 2);
            ComplexNumber c2 = new ComplexNumber(3, 4);

            // Act
            ComplexNumber sum = c1 + c2;

            // Assert
            Assert.AreEqual("4 + 6i", sum.Real + " + " + sum.Imaginary + "i");
        }

        [Test]
        public void ComplexNumber_Subtraction_Returns_CorrectResult()
        {
            // Arrange
            ComplexNumber c1 = new ComplexNumber(1, 2);
            ComplexNumber c2 = new ComplexNumber(3, 4);

            // Act
            ComplexNumber difference = c1 - c2;

            // Assert
            Assert.AreEqual("-2 + -2i", difference.Real.ToString() + " + " + difference.Imaginary.ToString() + "i");
        }

        [Test]
        public void ComplexNumber_Multiplication_Returns_CorrectResult()
        {
            // Arrange
            ComplexNumber c1 = new ComplexNumber(1, 2);
            ComplexNumber c2 = new ComplexNumber(3, 4);

            // Act
            ComplexNumber product = c1 * c2;

            // Assert
            Assert.AreEqual("-5 + 10i", product.Real.ToString() + " + " + product.Imaginary.ToString() + "i");
        }

        [Test]
        public void ComplexNumber_Division_Returns_CorrectResult()
        {
            // Arrange
            ComplexNumber c1 = new ComplexNumber(1, 2);
            ComplexNumber c2 = new ComplexNumber(3, 4);

            // Act
            ComplexNumber quotient = c1 / c2;

            // Assert
            Assert.AreEqual("0.44 + 0.08i", quotient.Real.ToString("0.00", CultureInfo.InvariantCulture) + " + " + quotient.Imaginary.ToString("0.00", CultureInfo.InvariantCulture) + "i");
        }

        [Test]
        public void ComplexNumber_Magnitude_Returns_CorrectMagnitude()
        {
            // Arrange
            ComplexNumber c1 = new ComplexNumber(3, 4);

            // Act
            double magnitude = c1.Magnitude();

            // Assert
            Assert.AreEqual(5, magnitude);
        }

        [Test]
        public void ComplexNumber_Phase_Returns_CorrectPhase()
        {
            // Arrange
            ComplexNumber c1 = new ComplexNumber(1, 1);

            // Act
            double phase = c1.Phase();

            // Assert
            Assert.AreEqual(Math.PI / 4, phase, 0.000001);
        }
    }
}