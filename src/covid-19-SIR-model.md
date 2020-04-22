# Modeling the Transmission Rate of COVID-19 with SIR in Smalltalk

## Using Glamorous Toolkit, PolyMath, Roassal2 and Ordinary Differential Equations

From [Wikipedia](https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology#The_SIR_model):
> The SIR model is one of the simplest compartmental models, and many models are derivatives of this basic form. The model consists of three compartments: S for the number of susceptible, I for the number of infectious, and R for the number of recovered or deceased (or immune) individuals. This model is reasonably predictive[citation needed] for infectious diseases that are transmitted from human to human, and where recovery confers lasting resistance, such as measles, mumps and rubella.

> These variables (S, I, and R) represent the number of people in each compartment at a particular time. To represent that the number of susceptible, infected and recovered individuals may vary over time (even if the total population size remains constant), we make the precise numbers a function of t (time): S(t), I(t) and R(t). For a specific disease in a specific population, these functions may be worked out in order to predict possible outbreaks and bring them under control. 

Lucky for us the [PolyMath](https://github.com/PolyMathOrg/PolyMath/) library contains everything we need to start solving Ordinary Differential Equations and even includes a [Tutorial](https://github.com/PolyMathOrg/PolyMath/wiki/Quick-start-to-ODE) on working with the SIR model.

From the wiki:

```smalltalk
|solver state system dt beta gamma values stepper diag|

dt := 1.0.
beta := 0.01.
gamma := 0.1.

system := PMExplicitSystem block: [ :x :t| |c|
	 c := Array new: 3.
	 c at: 1 put: (beta negated) * (x at: 1) * (x at: 2).
	 c at: 2 put: (beta * (x at: 1) * (x at: 2)) - (gamma * (x at: 2)).
	 c at: 3 put: gamma * (x at: 2).
	 c
	 ].

stepper := PMRungeKuttaStepper onSystem: system.
solver := (PMExplicitSolver new) stepper: stepper; system: system; dt: dt.
state := #(99 1 0).
values := (0.0 to: 200.0 by: dt) collect: [ :t| state := stepper doStep: state 
                                                          time: t stepSize: dt ].
```

In order to begin using this code to create our own model we must first understand the `dt` , `beta` and `gamma` variables.

- `dt` is a measure of time and in our model we'll use `1.0` so that each step represents `1 day`. 
- `beta` is the **transmission rate**. We'll discuss this further in the following section.
- `gamma` is the **recovery rate**. We'll discuss this further in the following section.

### The SIR Model

I found the best explanations of the SIR model to be this [Wikipedia](https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology#The_SIR_model) entry and the documentation of the R library [ShinySIR](https://cran.r-project.org/web/packages/shinySIR/vignettes/Vignette.html).

From `ShinySIR`:

> In the simple SIR model (without births or deaths), susceptible individuals (S) become infected and move into the infected class (I). After some period of time, infected individuals recover and move into the recovered (or immune) class (R). Once immune, they remain so for life (i.e. they do not leave the recovered class). The corresponding equations are given by:

>- \\( \frac{dS}{dt} = - \beta {S}{I} \\)
>- \\( \frac{dI}{dt} = \beta {S}{I} - \gamma I \\)
>- \\( \frac{dR}{dt} =  \gamma I \\)

> where `S`,`I`, and `R`, are the numbers of susceptible, infected, and recovered individuals in the population. Suppose the unit of time we are considering is days, then

>- β is the transmission rate and βSI represents the number of susceptible individuals that become infected per day;
>- γ is the recovery rate and γI is the number of infected individuals that recover per day;
>- 1/γ is the infectious period i.e. the average duration of time an individual remains infected.

> An important quantity of any disease model is the the reproductive number, R0, which represents the average number of secondary infections generated from one infectious individual in a completely susceptible population. For the SIR model,

> \\( R_0 = \beta N / \gamma, \\)

> where \\( N = S + I + R \\) is the total (constant) population size. Since R0 and the infectious period are more intuitive parameters, we use these as inputs for the built-in SIR model. We can then calculate β as:

> \\( \beta = R_0 \gamma / N. \\)

It's important to note that:

- \\( \beta \\) is the mathematical notation for `beta`.
- \\( \gamma \\) is the mathematical notation for `gamma`.

In order to calculate the `beta` parameter I source \\( R_0 \\) from [Report 3: Transmissibility of 2019-nCoV](https://www.imperial.ac.uk/media/imperial-college/medicine/sph/ide/gida-fellowships/Imperial-College-COVID19-transmissibility-25-01-2020.pdf) as `2.6`. This means that on average each infected person will infect 2.6 others. According to [CoronaTracker: World-wide COVID-19 Outbreak Data Analysis and Prediction](https://www.who.int/bulletin/online_first/20-255695.pdf) the **average duration of recovery** is `7.5 days`.

Given the average duration of recovery it's trivial to calculate our `gamma` value:

\\( \gamma = \frac{1}{D} \\) 

where \\( D \\) is the **average duration of recovery**.

```smalltalk
| durationOfRecovery gamma |

durationOfRecovery := 7.5.
gamma := 1 / durationOfRecovery
```

Given  \\( R_0 \\), \\( \gamma \\) and \\( N \\) (the total population size) it's trivial to calculate our `beta` value:

We'll start by setting `N` to the total population size of Florida. According to [Google](https://www.google.com/search?rlz=1C1CHBF_enUS880US880&sxsrf=ALeKk03BpRO2jyob9fJfErotMEWyBndmDw%3A1586996207554&ei=76OXXpmtIfKmggfIrJngDw&q=population+of+florida&oq=population+o&gs_lcp=CgZwc3ktYWIQARgAMgQIIxAnMgIIADICCAAyAggAMgIIADIFCAAQgwEyAggAMgIIADICCAAyAggAOgQIABBHOgUIABCRAjoECAAQQ0owCBcSLDBnNDA1ZzMzMGcyMDZnMTY2ZzE4MmcxNTJnMTU2ZzE0OWcxODNnMTgzZzE4ShsIGBIXMGcxZzFnMWcxZzFnMWcxZzFnMWczZzRQ1Pc1WN-DNmDrjDZoAHAEeACAAc4CiAGfEJIBCDAuMTAuMC4ymAEAoAEBqgEHZ3dzLXdpeg&sclient=psy-ab) the population of Florida is `21.48 million`.

```smalltalk
| durationOfRecovery gamma populationSize r0 beta |

durationOfRecovery := 7.5.
gamma := 1 / durationOfRecovery.
populationSize := 21480000.
r0 := 2.6.

beta := r0 * gamma / populationSize
```

We now have enough information to start simulating COVID-19 transmission across a population using SIR.

### Simulating COVID-19 transmission in Florida

Let's plug the information we have into the ODE tutorial from PolyMath.

```smalltalk
| dt durationOfRecovery gamma populationSize r0 beta solver state system values stepper |

dt := 1.0.
durationOfRecovery := 7.5.
gamma := 1 / durationOfRecovery.
populationSize := 21480000.
r0 := 2.6.
beta := r0 * gamma / populationSize.

system := PMExplicitSystem block: [ :x :t| |c|
	 c := Array new: 3.
	 c at: 1 put: (beta negated) * (x at: 1) * (x at: 2).
	 c at: 2 put: (beta * (x at: 1) * (x at: 2)) - (gamma * (x at: 2)).
	 c at: 3 put: gamma * (x at: 2).
	 c
	 ].

stepper := PMRungeKuttaStepper onSystem: system.
solver := (PMExplicitSolver new) stepper: stepper; system: system; dt: dt.
state := #(21480000 1 0).
values := (0.0 to: 180.0 by: dt) collect: [ :t| state := stepper doStep: state 
                                                          time: t stepSize: dt ].
```

This results in a multi-dimensional array indexed by days from the start of the spread to the 180th day (~6 months.) Each entry in the array is an array containing the Susceptible, Infected, and Recovered individuals. This information is much easier for humans to digest when visualized. For this we'll use the wonderful [Roassal2](http://agilevisualization.com/).

Let's take the values from our... uhhhh... *values* array and build some Roassal Data Frames.

```smalltalk
| dt durationOfRecovery gamma populationSize r0 beta solver state system values stepper susceptible infected recovered b ds1 ds2 ds3 |

dt := 1.0.
durationOfRecovery := 7.5.
gamma := 1 / durationOfRecovery.
populationSize := 21480000.
r0 := 2.6.
beta := r0 * gamma / populationSize.

system := PMExplicitSystem block: [ :x :t| |c|
	 c := Array new: 3.
	 c at: 1 put: (beta negated) * (x at: 1) * (x at: 2).
	 c at: 2 put: (beta * (x at: 1) * (x at: 2)) - (gamma * (x at: 2)).
	 c at: 3 put: gamma * (x at: 2).
	 c
	 ].

stepper := PMRungeKuttaStepper onSystem: system.
solver := (PMExplicitSolver new) stepper: stepper; system: system; dt: dt.
state := #(21480000 1 0).
values := (0.0 to: 180.0 by: dt) collect: [ :t| state := stepper doStep: state 
                                                          time: t stepSize: dt ].

susceptible := values collectWithIndex: [ :each :idx | Point x: idx y: (each at: 1) ].
infected := values collectWithIndex: [ :each :idx | Point x: idx y: (each at: 2) ].
recovered := values collectWithIndex: [ :each :idx | Point x: idx y: (each at: 3) ].

b := RTGrapher new.

ds1 := RTData new.
ds1 label: 'Susceptible'.
ds1 noDot.
ds1 points: susceptible.
ds1 connectColor: Color blue.
ds1 y: [ :v | v y ].
ds1 x: [ :v | v x ].
ds1 interaction popup text: [ :v | ((v key y) asInteger) asString, ' susceptible on day ', (v key x asString)].
b add: ds1.

ds2 := RTData new.
ds2 label: 'Infected'.
ds2 noDot.
ds2 points: infected.
ds2 connectColor: Color red.
ds2 y: [ :v | v y ].
ds2 x: [ :v | v x ].
ds2 interaction popup text: [ :v | ((v key y) asInteger) asString, ' infected on day ', (v key x asString)].
b add: ds2.

ds3 := RTData new.
ds3 label: 'Recovered'.
ds3 noDot.
ds3 points: recovered.
ds3 connectColor: Color green.

ds3 y: [ :v | v y ].
ds3 x: [ :v | v x ].
ds3  interaction popup text: [ :v | ((v key y) asInteger) asString, ' recovered on day ', (v key x asString)].
b add: ds3.

b addDecorator: (RTCursorFollower new color: Color gray).
b axisX title: 'Days'; noDecimal.
b axisY title: 'Population'; noDecimal.
b legend right.

b build view canvas buildMorph extent: 1000@500; exportAsPNG.
```

![Florida COVID-19 transmission](https://i.imgur.com/6V1lpix.png)

As you can see, this lays out our data in an easily digestable fashion. It shows the peak of the infections happening on **day 84** with **5,143,437 individuals infected**.