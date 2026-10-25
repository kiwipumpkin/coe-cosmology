Residuals as Features: A Framework for Scientific Discovery

1. Central Thesis
Residuals from a sufficiently accurate model should be treated as observable features worthy of scientific investigation, rather than merely as errors to minimize.
Every scientific model begins with assumptions, predictions, simplifications, and a limited selection of variables.
To the extent that the model accurately describes the phenomena it incorporates, the remaining residuals represent aspects of the observations that the model does not explain.
These may include phenomena that:
- Were deliberately excluded for simplicity.
- Are known but imperfectly modeled.
- Were not recognized as relevant.
- Represent previously unidentified physical mechanisms.
- Arise from measurement error or irreducible stochastic behavior.
The central proposition is that the unexplained remainder of a successful model can itself become an object of scientific investigation.

2. Mathematical Formulation
Consider an observed system: \[ Y=M(X,\theta)+R \]
where: Y: Observations. M(X,θ): Predictions from the existing model. R: Residuals representing unexplained variation.
The residuals are: \[ \boxed{R=Y-M(X,\theta)} \]
Traditional model optimization often seeks to minimize: \[ \min_\theta \sum_i R_i^2 \]
However, scientific discovery introduces a complementary objective: \[ \boxed{\text{Identify and characterize structure within }R} \]
These objectives are related but distinct.
A model can produce small residuals while leaving scientifically important systematic patterns unexplained.
Conversely, a model with relatively large residuals may be scientifically useful if those residuals are understood and characterized.
Residual magnitude measures unexplained deviation; residual structure may provide information about its origin.

3. Residuals as Observable Features
Residuals should be examined for measurable properties just as the original observations are.
TABLE (Residual feature | Potential significance):
Periodicity | Missing oscillatory processes
Spatial correlation | Unmodeled spatial interactions
Directional dependence | Anisotropy or broken symmetry
Temporal correlation | Missing dynamics or memory effects
Systematic curvature | Incorrect functional relationship
Discontinuities | Regime changes or phase transitions
Changing variance | Heteroscedasticity or evolving processes
Spectral structure | Missing characteristic scales or frequencies
Statistical independence | Consistency with noise under the tested assumptions
A residual is not automatically evidence of new physics.
Nevertheless, reproducible structure that survives appropriate statistical testing may indicate missing explanatory variables or mechanisms.
The appropriate question becomes:
What physical, mathematical, or observational process could generate the structure remaining after the known contributions have been modeled?

4. Progressive Isolation of Unexplained Phenomena
Suppose a physical system contains several contributions: \[ Y=F_1+F_2+F_3+F_4+\epsilon \]
An initial model incorporates only two: \[ M_1=F_1+F_2 \]
Its residuals are: \[ R_1=F_3+F_4+\epsilon \]
After identifying and incorporating another contribution: \[ M_2=F_1+F_2+F_3 \]
The residual becomes: \[ R_2=F_4+\epsilon \]
As established phenomena are accounted for, previously obscured contributions may become distinguishable.
This suggests an iterative discovery process: \[ \boxed{\text{Observation}\rightarrow\text{Model}\rightarrow\text{Residual analysis}\rightarrow\text{Hypothesis}\rightarrow\text{Revised model}} \]
The process can repeat as new physical mechanisms are identified.
Importantly, this decomposition is an idealization. Physical contributions may interact nonlinearly, and statistical fitting may absorb or distort signals from omitted variables.

5. Historical Examples
Several important scientific discoveries illustrate this approach.
Johannes Kepler — Planetary Orbits
Kepler encountered an approximately eight-arcminute discrepancy between observations of Mars and the circular orbital models he investigated.
Rather than dismissing the discrepancy, he pursued an alternative geometric description.
This contributed to the development of elliptical planetary orbits.
Lesson: Residuals can reveal incorrect assumptions about the mathematical structure of a physical system.
Urbain Le Verrier — Discovery of Neptune
Unexplained deviations in Uranus's orbit motivated calculations predicting an additional planet.
Neptune was subsequently observed near its predicted location.
Lesson: Residuals can reveal missing physical objects without requiring replacement of the underlying governing theory.
Albert Einstein — Mercury's Perihelion
After accounting for known Newtonian gravitational effects, approximately 43 arcseconds per century of Mercury's perihelion precession remained unexplained.
General relativity subsequently accounted for this discrepancy.
Lesson: Residuals can reveal limitations in otherwise highly successful physical theories.
These examples demonstrate that residual analysis can lead to different classes of discovery:
- Previously unrecognized objects.
- Previously omitted interactions.
- Incorrect mathematical assumptions.
- New physical laws.

6. Researchers and Intellectual Foundations
TABLE (Researcher | Contribution | Relevance):
John Tukey | Exploratory data analysis | Residuals as sources of discoverable structure
George E. P. Box | Iterative scientific modeling and model criticism | Learning through systematic model inadequacies
Aris Spanos | Statistical misspecification testing | Identifying violated modeling assumptions
Deborah Mayo | Severe testing and scientific inference | Distinguishing genuine discrepancies from statistical artifacts
Johannes Kepler | Planetary orbital analysis | Residuals motivating new mathematical descriptions
Urbain Le Verrier | Prediction of Neptune | Residuals indicating omitted physical entities
Albert Einstein | General relativity | Residuals revealing limitations of established physical laws
Particularly Relevant Literature
John Tukey — Exploratory Data Analysis (1977): Emphasizes discovering structure through graphical and statistical examination of observations, including residuals.
George Box — Science and Statistics (1976): Develops the view of scientific modeling as an iterative process involving construction, criticism, and revision.
Aris Spanos — Work on Statistical Misspecification: Investigates how systematic departures from model assumptions can be identified and evaluated.
Deborah Mayo — Error and the Growth of Experimental Knowledge (1996): Examines how rigorous testing of discrepancies supports scientific inference.
These works provide complementary foundations for treating residual analysis as part of scientific discovery.

7. Distinguishing Prediction from Discovery
Two different objectives emerge.
Predictive modeling: The primary objective is to minimize disagreement between observations and predictions. \[ \min_\theta L(Y,M) \] where L is an appropriate loss function.
A sufficiently flexible model may achieve excellent predictive accuracy without providing a useful physical explanation.
Discovery-oriented modeling: The objective is to explain established phenomena while preserving the ability to recognize previously unidentified structure.
A simplified conceptual formulation is: \[ \boxed{\text{Scientific discovery}=\text{Accurate modeling}+\text{Systematic investigation of residuals}} \]
This does not imply that large residuals are desirable.
Rather, it recognizes that minimizing residual magnitude and understanding residual structure are different scientific activities.
Overfitting is particularly problematic because it can remove genuine signals before they are recognized as distinct phenomena.

8. Philosophical Implications
The discussion suggests several broader principles.
First: A model need not be complete to be scientifically successful. A model can accurately represent the phenomena it was designed to describe while intentionally excluding other contributions.
Second: Residuals are defined relative to the model. They do not constitute an objective separation between known and unknown reality. Different models can produce different residual structures from the same observations.
Third: The unexplained remainder can have independent scientific value. Once calculated, residuals can be investigated for correlations, symmetries, characteristic scales, distributions, and temporal evolution.
Fourth: A sufficiently accurate model can function as a filter. By accounting for established effects, it may expose phenomena that were previously obscured.
Fifth: Discovery requires distinguishing unexplained structure from observational artifacts. A reproducible residual pattern is a candidate for explanation, not automatically evidence of a new physical mechanism.

9. Proposed Research Framework
The discussion suggests a systematic methodology:
1. Construct a defensible baseline model. Explicitly identify assumptions, simplifications, and omitted effects.
2. Validate the model. Establish its predictive accuracy and characterize measurement uncertainty.
3. Calculate residuals. Preserve their relationships to relevant observational variables.
4. Analyze residual features. Investigate spatial, temporal, spectral, directional, and statistical structure.
5. Test against known omissions. Determine whether established but excluded mechanisms explain the observed features.
6. Develop alternative hypotheses. Propose physical mechanisms capable of generating the remaining structure.
7. Test independently. Evaluate the hypotheses against new observations, withheld data, or independent measurements.
8. Revise the model. Incorporate supported mechanisms and repeat the analysis.
This framework combines predictive modeling, exploratory analysis, statistical testing, and hypothesis generation.

10. Principal Conclusion
The discussion began with the observation that every model necessarily incorporates assumptions and simplifications.
As a model becomes more accurate, its residuals increasingly characterize the phenomena it does not explain, although they may also contain noise, measurement errors, and artifacts of the modeling procedure.
Consequently, residuals should not be regarded exclusively as undesirable deviations to eliminate.
They can be treated as measurable features, investigated systematically, and used to motivate new scientific hypotheses.
The distinction is between using residuals to judge the quality of a model and using residuals to investigate phenomena beyond the model's explanatory scope.
Consolidated Thesis (block quote):
Every scientific model is necessarily selective in the phenomena it represents. To the extent that a model accurately accounts for its intended physical contributions, its residuals characterize the unexplained remainder of the observations. Rather than treating these residuals exclusively as errors to minimize, researchers can investigate their structure as observable features. Systematic, reproducible patterns may reveal omitted variables, unrecognized interactions, incorrect assumptions, or previously unidentified phenomena. In this sense, a successful model serves not only as an explanation of what is understood, but also as an instrument for investigating what remains unexplained.
Central proposition: The residuals of a good model are not merely measures of its limitations; they are potential observations of phenomena beyond its current explanatory scope.
