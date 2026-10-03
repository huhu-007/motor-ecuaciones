"""Suite de tests del motor (unittest, sin dependencias).  Ejecutar:  python -m unittest -v tests.test_engine"""
import time
import unittest
from fractions import Fraction as F

from mpmath import mp, mpf

from equation_engine import EquationEngine, parse_equation, ParseError
from equation_engine.parser import ExpressionParser
from equation_engine.nodes import to_text
from equation_engine.evaluator import FunctionEvaluator
from equation_engine.backends import MPBackend, IntervalBackend
from equation_engine.exact import ExactBackend, Rad
from equation_engine.errors import DomainError, Uncertified, NotExact
from equation_engine.domain import DomainAnalyzer
from equation_engine.bisection import buscar_periodo, fraccion_periodica
from equation_engine.polynomial import reduced_ratfunc, yun, poly_text
from equation_engine.bolzano import BOLZANO_NO_APLICABLE

ENG = EquationEngine(digits=20, analysis_digits=50, analysis_budget=6.0)


def fx(text):
    return FunctionEvaluator(parse_equation(text).function)


def mpval(text, x):
    be = MPBackend(40)
    return fx(text).evaluate(be, be.from_fraction(F(x)))


def solve(t):
    return ENG.solve(t)


def values(r):
    return sorted(float(x["decimal_form"]) for x in r["roots"])


class ParserTests(unittest.TestCase):
    VALID = ["x", "2x", "3(x+1)", "(x+1)(x-1)", "2π", "2sin(x)", "x^2", "x^3", "x^(x+1)", "2^x", "e^x", "√x",
             "√(x+1)", "∛x", "⁴√x", "¹⁰⁰√x", "root_7(x)", "log(x)", "log_2(x)", "log_10(x)", "ln(x)",
             "ln(x+1)", "log_3(x^2)", "π/2", "e^(2x)", "x^x", "1/(√2)", "√3/√5", "(x+1)/(x-2)",
             "sin(x)+x", "2cos(x)-1", "e^x-sin(x)", "1,5x", "x²", "sin^2(x)", "tan(x)", "x^3 = 2x+5", "-x^2", "2^-x"]

    def test_valid_expressions_parse(self):
        for t in self.VALID:
            with self.subTest(t=t):
                parse_equation(t)

    def test_text_roundtrip_same_value(self):
        for t in ["x^2-2", "2x(x+1)", "√(x+1)*2", "e^(2x)-sin(x)", "log_2(x^2+1)", "-x^2+3", "2^-x", "1/(x+2)"]:
            with self.subTest(t=t):
                a = parse_equation(t).function
                b = ExpressionParser().parse(to_text(a))
                x = F(13, 10)
                be = MPBackend(30)
                self.assertEqual(FunctionEvaluator(a).evaluate(be, be.from_fraction(x)),
                                 FunctionEvaluator(b).evaluate(be, be.from_fraction(x)))

    def test_precedence(self):
        self.assertEqual(mpval("-x^2", 3), -9)
        self.assertEqual(mpval("2^3^2", 0), 512)
        self.assertEqual(mpval("1/2x", 4), 2)           # (1/2)·x, como el código original
        self.assertEqual(mpval("2^-2", 0), mpf(1) / 4)

    def test_pi_e_are_constants_not_variables(self):
        self.assertEqual(parse_equation("e^x").variable, "x")
        self.assertEqual(parse_equation("2π").variable, "x")
        self.assertAlmostEqual(float(mpval("e*x", 1)), 2.718281828, places=8)       # e*x  != e^x
        self.assertAlmostEqual(float(mpval("e^x", 1)), 2.718281828, places=8)
        self.assertAlmostEqual(float(mpval("e^x", 2)), 7.389056099, places=8)

    def test_rejected_inputs(self):
        for t in ["x+", "", "sin 2x", "√x^2", "1 000", "3e5", "x=1=2", "sec(x)", "root_1(x)", "root_101(x)",
                  "x|y", "2 3", "x+y", "((x)", "x^", "@x"]:
            with self.subTest(t=t), self.assertRaises(ParseError):
                parse_equation(t)

    def test_security_no_code_execution(self):
        for t in ["__import__('os').system('echo hacked')", "open('f')", "x.__class__", "lambda: 1", "exec('1')",
                  "import os", "x; y"]:
            with self.subTest(t=t), self.assertRaises(ParseError):
                parse_equation(t)

    def test_too_long_and_deep(self):
        with self.assertRaises(ParseError):
            parse_equation("x+" * 400 + "1")
        with self.assertRaises(ParseError):
            parse_equation("(" * 200 + "x" + ")" * 200)


class EvaluatorTests(unittest.TestCase):
    def ex(self, t, x=0):
        return fx(t).evaluate(ExactBackend(), Rad.rat(x))

    def test_exact_roots(self):
        self.assertEqual(self.ex("√4").rational(), 2)
        self.assertEqual(self.ex("∛(-8)").rational(), -2)
        self.assertEqual(self.ex("⁴√16").rational(), 2)
        self.assertEqual(self.ex("¹⁰⁰√(2^100)").rational(), 2)

    def test_even_root_negative_undefined_odd_ok(self):
        with self.assertRaises(DomainError):
            self.ex("√(-4)")
        with self.assertRaises(DomainError):
            self.ex("⁴√(-16)")
        self.assertEqual(self.ex("root_5(-32)").rational(), -2)
        self.assertEqual(self.ex("root_99(-1)").rational(), -1)

    def test_log_domain(self):
        for t in ["ln(0)", "ln(-1)", "log(0)", "log_2(-4)", "log_1(5)", "log_(-2)(4)", "log_0(3)"]:
            with self.subTest(t=t), self.assertRaises(DomainError):
                be = MPBackend(30)
                fx(t).evaluate(be, be.from_fraction(F(0)))
        self.assertEqual(self.ex("log_2(8)").rational(), 3)
        self.assertEqual(self.ex("log_10(1000)").rational(), 3)

    def test_division_by_zero(self):
        with self.assertRaises(DomainError):
            self.ex("1/(x-2)", 2)
        with self.assertRaises(DomainError):
            mpval("(x+1)/(x-2)", 2)

    def test_tan_undefined_at_pi_half(self):
        be = ExactBackend()
        half_pi = be.div(be.const("pi"), Rad.rat(2))
        with self.assertRaises(DomainError):
            fx("tan(x)").evaluate(be, half_pi)
        self.assertEqual(fx("tan(x)").evaluate(be, Rad.rat(0)).rational(), 0)

    def test_power_conventions(self):
        self.assertEqual(self.ex("x^(1/3)", -8).rational(), -2)         # raíz impar real
        with self.assertRaises(DomainError):
            self.ex("x^(1/2)", -4)
        with self.assertRaises(DomainError):
            self.ex("x^x", -1)
        with self.assertRaises(DomainError):
            self.ex("0^0")
        self.assertEqual(self.ex("2^x", 3).rational(), 8)

    def test_exact_radical_arithmetic(self):
        r = self.ex("√2/√3")
        self.assertEqual(r.t, Rad._factor_fraction_root(F(2, 3), 2).t)
        z = self.ex("(√3+√2)^2-5-2√6")
        self.assertTrue(z.is_zero())
        self.assertTrue(self.ex("(1+√2)(1-√2)+1").is_zero())

    def test_interval_rejects_poles(self):
        iv = IntervalBackend(30)
        f = fx("1/(x-2)")
        with self.assertRaises(DomainError):
            f.evaluate(iv, iv.interval(F(19, 10), F(21, 10)))
        f.evaluate(iv, iv.interval(F(1), F(19, 10)))                    # sin polo: certificable
        t = fx("tan(x)")
        with self.assertRaises(DomainError):
            t.evaluate(iv, iv.interval(F(3, 2), F(8, 5)))               # contiene π/2


class DomainTests(unittest.TestCase):
    def dom(self, t):
        d = DomainAnalyzer()
        return d.text(d.analyze(parse_equation(t).function))

    def test_domains(self):
        self.assertEqual(self.dom("1/(x-2)"), "x ≠ 2")
        self.assertEqual(self.dom("ln(x)"), "x > 0")
        self.assertEqual(self.dom("log(x)"), "x > 0")
        self.assertEqual(self.dom("√(x-3)"), "x ≥ 3")
        self.assertIn("π/2", self.dom("tan(x)"))
        self.assertEqual(self.dom("log_3(x^2)"), "x ≠ 0")
        self.assertEqual(self.dom("x^x"), "x > 0")
        self.assertEqual(self.dom("2^x"), "ℝ (sin restricciones)")
        self.assertEqual(self.dom("∛x"), "ℝ (sin restricciones)")
        self.assertEqual(self.dom("⁴√(2x-1)"), "x ≥ 1/2")
        self.assertIn("x ≠ 0", self.dom("log_(x+1)(x)"))


class EquationTests(unittest.TestCase):
    """Los casos pedidos en el enunciado."""

    def test_x2_minus_2(self):
        r = solve("x^2 - 2 = 0")
        self.assertEqual(r["status"], "solved")
        self.assertEqual([x["exact_form"] for x in r["roots"]], ["-√2", "√2"])
        self.assertTrue(all(x["exactness"] == "exact_proven" for x in r["roots"]))
        self.assertAlmostEqual(values(r)[1], 2 ** 0.5, places=12)
        self.assertTrue(r["bolzano_applicable"])

    def test_x3_minus_8(self):
        r = solve("x^3 - 8 = 0")
        self.assertEqual(r["roots"][0]["exact_form"], "2")
        self.assertEqual(r["roots"][0]["exactness"], "exact_proven")
        self.assertEqual(len(r["roots"]), 1)

    def test_exp(self):
        r = solve("e^x - 2 = 0")
        self.assertEqual(r["status"], "solved")
        x = r["roots"][0]
        self.assertAlmostEqual(float(x["decimal_form"]), 0.6931471805599453, places=14)
        self.assertEqual(x["candidate_form"], "ln(2)")
        self.assertEqual(x["exactness"], "numerical_evidence")        # ln(2) no se demuestra con el backend exacto
        self.assertIsNone(x["exact_form"])                            # y NO se anuncia como exacta

    def test_ln(self):
        r = solve("ln(x) - 1 = 0")
        x = r["roots"][0]
        self.assertAlmostEqual(float(x["decimal_form"]), 2.718281828459045, places=13)
        self.assertEqual(x["exact_form"], "e")
        self.assertEqual(x["exactness"], "exact_verified")
        self.assertEqual(r["domain"]["text"], "x > 0")

    def test_sin(self):
        r = solve("sin(x) = 0")
        v = values(r)
        self.assertEqual(len(v), 63)                                   # kπ, k = -31..31
        self.assertTrue(any(abs(a) < 1e-12 for a in v))
        self.assertTrue(any(abs(a - 3.141592653589793) < 1e-12 for a in v))
        near = [x for x in r["roots"] if abs(float(x["decimal_form"]) - 3.141592653589793) < 1e-9][0]
        self.assertEqual(near["exact_form"], "π")

    def test_cos_minus_1_even_multiplicity_roots(self):
        """cos(x)-1 ≤ 0: NO hay cambio de signo; Bolzano no ve las raíces 2πk, pero se verifican exactamente."""
        r = solve("cos(x) - 1 = 0")
        self.assertEqual(r["status"], "solved")
        forms = {x["exact_form"] for x in r["roots"]}
        self.assertIn("0", forms)
        self.assertIn("2π", forms)
        self.assertIn("-2π", forms)
        tang = [x for x in r["roots"] if x["method"] == "tangent_root_exact"]
        self.assertTrue(tang)
        self.assertFalse(tang[0]["bolzano"]["applicable"])

    def test_linear(self):
        r = solve("x/2 - 3 = 0")
        self.assertEqual(r["roots"][0]["exact_form"], "6")
        self.assertEqual(r["roots"][0]["exactness"], "exact_proven")

    def test_sqrt(self):
        r = solve("√x - 2 = 0")
        self.assertEqual(r["roots"][0]["exact_form"], "4")
        self.assertEqual(r["domain"]["text"], "x ≥ 0")

    def test_cbrt(self):
        r = solve("∛x + 2 = 0")
        self.assertEqual(r["roots"][0]["exact_form"], "-8")

    def test_root_found_near_domain_boundary(self):
        r = solve("ln(x) + 5 = 0")
        self.assertAlmostEqual(float(r["roots"][0]["decimal_form"]), 0.006737946999085467, places=15)
        self.assertEqual(r["roots"][0]["exact_form"], "e^-5")

    def test_rational_via_periodicity_then_exact_proof(self):
        r = solve("3x - 1 = 0")
        x = r["roots"][0]
        self.assertEqual(x["exact_form"], "1/3")
        self.assertEqual(x["exactness"], "exact_proven")
        self.assertIn("(3)", x["periodicity"]["period_notation"])


class BolzanoTests(unittest.TestCase):
    def test_double_root_x2_is_found(self):
        r = solve("x^2 = 0")
        self.assertEqual(r["status"], "solved")
        self.assertEqual(r["roots"][0]["exact_form"], "0")
        self.assertEqual(r["roots"][0]["multiplicity"], 2)

    def test_double_root_not_on_grid(self):
        r = solve("(x - 1/3)^2 = 0")
        self.assertEqual(r["roots"][0]["exact_form"], "1/3")
        self.assertEqual(r["roots"][0]["multiplicity"], 2)

    def test_tangent_root_with_radical_non_polynomial(self):
        r = solve("(x - √2)^2 = 0")
        self.assertEqual(r["status"], "solved")
        self.assertEqual(r["roots"][0]["exact_form"], "√2")
        self.assertEqual(r["roots"][0]["method"], "tangent_root_exact")

    def test_poles_are_not_roots(self):
        r = solve("tan(x) = 0")
        self.assertEqual(len(r["roots"]), 63)
        self.assertGreaterEqual(len(r["rejected_sign_changes"]), 60)
        poles = [x for x in r["roots"] if abs(abs(float(x["decimal_form"]) % 3.141592653589793) - 1.5707963) < 0.01]
        self.assertEqual(poles, [])

    def test_1_over_x_has_no_solution(self):
        r = solve("1/x = 0")
        self.assertEqual(r["status"], "no_real_solution")
        self.assertEqual(r["roots"], [])
        self.assertTrue(r["proven_no_real_solution"])

    def test_1_over_x_generic_route_rejects_pole(self):
        r = solve("1/x + sin(0)*x = 0")                  # fuerza la ruta general (sin racional)
        self.assertEqual(r["roots"], [])
        self.assertNotEqual(r["status"], "solved")

    def test_removable_discontinuity_excluded(self):
        r = solve("(x^2 - 1)/(x - 1) = 0")
        self.assertEqual([x["exact_form"] for x in r["roots"]], ["-1"])
        self.assertEqual(r["domain"]["text"], "x ≠ 1")

    def test_pole_with_numerator_root(self):
        r = solve("(x + 1)/(x - 2) = 0")
        self.assertEqual([x["exact_form"] for x in r["roots"]], ["-1"])

    def test_no_real_solution_proven_polynomial(self):
        r = solve("x^2 + 1 = 0")
        self.assertEqual(r["status"], "no_real_solution")
        self.assertTrue(r["proven_no_real_solution"])

    def test_roots_outside_range_are_not_no_solution(self):
        r = solve("x^2 - 40000 = 0")
        self.assertEqual(r["status"], "no_root_found_in_range")
        self.assertFalse(r["proven_no_real_solution"])
        self.assertTrue(any("FUERA" in w for w in r["warnings"]))

    def test_undefined_in_range(self):
        r = solve("√(-1 - x^2) = 0")
        self.assertEqual(r["status"], "undefined_in_range")

    def test_bolzano_message_when_not_applicable(self):
        """Raíz tangente no identificable: se debe decir que Bolzano no se puede aplicar,
        NO que la ecuación no tenga solución."""
        r = solve("(x - sin(1)/3)^2 = 0")
        self.assertEqual(r["status"], "bolzano_not_applicable")
        self.assertEqual(r["bolzano_message"], BOLZANO_NO_APLICABLE)
        self.assertIn(BOLZANO_NO_APLICABLE, r["messages"])
        self.assertFalse(r["proven_no_real_solution"])
        self.assertFalse(r["bolzano_applicable"])
        self.assertTrue(r["unverified_candidates"])

    def test_proven_no_root_in_range_is_scoped(self):
        r = solve("sin(x)^2 + 1/1000 = 0")
        self.assertEqual(r["status"], "no_root_found_in_range")
        self.assertTrue(r["diagnostics"].get("proven_no_root_in_range"))
        self.assertFalse(r["proven_no_real_solution"])      # solo dentro del rango

    def test_identity_and_constants(self):
        self.assertEqual(solve("x = x")["status"], "identity")
        self.assertEqual(solve("1 = 2")["status"], "no_real_solution")
        self.assertEqual(solve("2 = 2")["status"], "identity")
        self.assertEqual(solve("√(-1) = 0")["status"], "undefined_in_range")

    def test_inequalities_parse_but_are_not_solved(self):
        r = solve("x^2 > 2")
        self.assertEqual(r["status"], "unsupported_relation")
        self.assertFalse(r["success"])

    def test_invalid_input_status(self):
        for t in ["x+", "sec(x)", "__import__('os')"]:
            r = solve(t)
            self.assertEqual(r["status"], "invalid_input", t)
            self.assertFalse(r["success"])

    def test_multiple_roots_polynomial_with_sturm_completeness(self):
        r = solve("(x-1)(x-2)(x-3)(x+1/2) = 0")
        self.assertEqual([x["exact_form"] for x in r["roots"]], ["-1/2", "1", "2", "3"])
        r = solve("x^2*(x-1/3)^2 = 0")
        self.assertEqual([x["exact_form"] for x in r["roots"]], ["0", "1/3"])

    def test_close_roots_in_one_grid_cell(self):
        r = solve("(x - 0.31)(x - 0.33) = 0")        # ambas dentro de la celda [0.3, 0.4]
        self.assertEqual([x["exact_form"] for x in r["roots"]], ["31/100", "33/100"])

    def test_irreducible_cubic_minimal_polynomial_proven(self):
        r = solve("x^3 - 2x - 5 = 0")
        x = r["roots"][0]
        self.assertAlmostEqual(float(x["decimal_form"]), 2.0945514815423265, places=14)
        self.assertEqual(x["exactness"], "none")
        self.assertEqual(x["minimal_polynomial"]["status"], "exact_proven")
        self.assertIn("x^3 - 2x - 5", x["minimal_polynomial"]["polynomial"])

    def test_high_exponent_does_not_break(self):
        r = solve("x^1000 - 1 = 0")
        self.assertEqual(sorted(x["exact_form"] for x in r["roots"]), ["-1", "1"])


def dec(v, n):
    """Cadena decimal con exactamente n decimales TRUNCADOS (como produce la bisección)."""
    t = mp.nstr(v, n + 15, strip_zeros=False)
    ip, _, fp = t.partition(".")
    return ip + "." + fp[:n]


class NumberAnalysisTests(unittest.TestCase):
    """√2, √3+√2, √2/√3, ⁴√16, ¹⁰⁰√(...)"""

    def test_sqrt2(self):
        r = ENG.analyze_expression("√2")
        self.assertEqual((r["exact_form"], r["exactness"]), ("√2", "exact_proven"))

    def test_sum_of_radicals(self):
        r = ENG.analyze_expression("√3 + √2")
        self.assertEqual((r["exact_form"], r["exactness"]), ("√2 + √3", "exact_proven"))

    def test_quotient(self):
        r = ENG.analyze_expression("√2 / √3")
        self.assertEqual(r["exact_form"], "√6/3")

    def test_fourth_root(self):
        r = ENG.analyze_expression("⁴√16")
        self.assertEqual(r["exact_form"], "2")

    def test_hundredth_root(self):
        r = ENG.analyze_expression("¹⁰⁰√(7/3)")
        self.assertEqual(r["exact_form"], "¹⁰⁰√(7/3)")
        r = ENG.analyze_expression("¹⁰⁰√2")
        self.assertEqual(r["exact_form"], "¹⁰⁰√2")
        r = ENG.analyze_expression("2*root_100(5)")
        self.assertEqual(r["exactness"], "exact_proven")

    def test_nested_and_products(self):
        self.assertEqual(ENG.analyze_expression("∛2 * √3")["exact_form"], "⁶√108")
        self.assertEqual(ENG.analyze_expression("(1+√2)^2")["exact_form"], "3 + 2√2")
        self.assertEqual(ENG.analyze_expression("1/(1+√2)")["exact_form"], "-1 + √2")

    def test_constants_kept_symbolic(self):
        self.assertEqual(ENG.analyze_expression("2π")["exact_form"], "2π")
        self.assertEqual(ENG.analyze_expression("π/2")["exact_form"], "π/2")
        self.assertEqual(ENG.analyze_expression("e^2*e")["exact_form"], "e^3")

    def test_not_exact_is_labelled_numerical(self):
        r = ENG.analyze_expression("sin(1) + √2")
        self.assertNotEqual(r["exactness"], "exact_proven")
        self.assertEqual(r["status"], "numerical_only")

    def test_even_root_of_negative_undefined(self):
        self.assertEqual(ENG.analyze_expression("√(-4)")["status"], "undefined")
        self.assertEqual(ENG.analyze_expression("∛(-8)")["exact_form"], "-2")

    # ---- reconocimiento NUMÉRICO (solo evidencia)
    def test_numeric_recognition_is_evidence_not_proof(self):
        mp.dps = 120
        rep = ENG.analyze_number(dec(mp.sqrt(2), 78), 78)
        self.assertEqual(rep["best"]["exact_form"], "√2")
        self.assertEqual(rep["best"]["level"], "numerical_evidence")
        self.assertIn("no es una demostración", rep["statement"])
        self.assertTrue(rep["square_is_rational"])

    def test_numeric_sum_nested_and_high_degree(self):
        mp.dps = 200
        s = lambda v, n: dec(v, n)
        self.assertEqual(ENG.analyze_number(s(mp.sqrt(3) + mp.sqrt(2), 78), 78)["best"]["exact_form"], "√2 + √3")
        self.assertEqual(ENG.analyze_number(s(mp.sqrt(1 + mp.sqrt(2)), 78), 78)["best"]["exact_form"], "√(1 + √2)")
        self.assertEqual(ENG.analyze_number(s(mp.root(2, 100), 118), 118)["best"]["exact_form"], "¹⁰⁰√2")
        self.assertEqual(ENG.analyze_number(s(mp.root(mp.mpf(7) / 3, 100), 118), 118)["best"]["exact_form"], "¹⁰⁰√(7/3)")
        self.assertEqual(ENG.analyze_number(s(mp.root(16, 4), 58), 58)["best"]["exact_form"], "2")

    def test_regression_digits_cannot_exceed_string(self):
        rep = ENG.analyze_number("2.0", 58)               # antes: 2.0 + 0.05 -> 41/20 (falso)
        self.assertNotEqual((rep["best"] or {}).get("exact_form"), "41/20")
        self.assertTrue(any("solo trae" in n for n in rep["notes"]))

    # ---- regresión de falsos positivos descubiertos durante el desarrollo
    def test_regression_no_huge_root_for_sqrt3_plus_sqrt2(self):
        mp.dps = 200
        rep = ENG.analyze_number(dec(mp.sqrt(3) + mp.sqrt(2), 80), 80)
        self.assertEqual(rep["best"]["exact_form"], "√2 + √3")
        self.assertLess(len(rep["best"]["exact_form"]), 20)        # nunca ²⁰√(número de 300 cifras)

    def test_regression_near_miss_is_not_identified(self):
        mp.dps = 200
        v = mp.sqrt(2) + mp.mpf(10) ** -12
        rep = ENG.analyze_number(dec(v, 40), 40)
        self.assertIsNone(rep["best"])                              # NO es √2
        self.assertFalse((rep["polynomial"] or {}).get("found"))    # ni un polinomio espurio de grado 7

    def test_regression_insufficient_precision_is_reported(self):
        mp.dps = 50
        rep = ENG.analyze_number(dec(mp.sqrt(2), 10), 10)
        self.assertTrue(rep["notes"])                               # avisa de precisión insuficiente
        # con tan pocos dígitos no se debe afirmar una combinación de radicales exótica
        self.assertTrue(rep["best"] is None or rep["best"]["level"] == "numerical_evidence")

    def test_regression_random_real_not_identified(self):
        mp.dps = 120
        v = mp.mpf("0.8472130847939790866064991234821223")          # número arbitrario
        rep = ENG.analyze_number(dec(v * 3, 48), 48)
        self.assertIsNone(rep["best"])


class OutputTests(unittest.TestCase):
    def test_json_serialisable_and_schema(self):
        import json
        from equation_engine.formatter import ResultFormatter
        r = solve("x^2 - 2 = 0")
        data = json.loads(ResultFormatter.to_json(r))
        for k in ["success", "status", "equation", "roots", "method", "bolzano_applicable", "domain", "warnings"]:
            self.assertIn(k, data)
        self.assertNotIn("_sortkey", json.dumps(data))
        self.assertIn("√2", ResultFormatter.to_text(r))

    def test_periodicity_helpers_original_behaviour(self):
        self.assertEqual(buscar_periodo([1, 2] + [3] * 30), ([1, 2], [3]))
        self.assertEqual(fraccion_periodica("", 0, ([], [3])), F(1, 3))
        self.assertEqual(fraccion_periodica("-", 0, ([1], [6])), F(-1, 6))
        self.assertIsNone(buscar_periodo([1, 4, 1, 4, 2, 1, 3, 5, 6, 2, 3, 7, 3, 0, 9, 5, 0, 4, 8, 8]))


if __name__ == "__main__":
    unittest.main(verbosity=2)
