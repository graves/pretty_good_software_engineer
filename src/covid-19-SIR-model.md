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

Our model is in no way precise and leaves out many factors such as the effects of social distancing, state measures to slow the spread, population density, age, etc. Even so, the peak infected rate seems quite high at over 5 million individuals while at the time of writing (midnight April 24, 2020) there are only 30,839 [confirmed infected individuals in Florida.](https://floridahealthcovid19.gov/).

Let's do a little investigating to see what the problem may be. We'll first change the Y-axis to our actual calendar date. The CDC [confirmed the first two COVID-19 cases in Florida](https://www.clickorlando.com/health/2020/03/02/gov-desantis-declares-health-emergency-after-2-presumptive-positive-coronavirus-cases-found/) on March 1, 2020 so we'll change our Y-axis to start on March 1. The linked article also announces 2 initial cases so we'll change our starting state so that it repressents 2 infected individuals.

`state := #(21480000 2 0).`

![Florida COVID-19 transmission infected](https://i.imgur.com/NY0aFUm.png)

The graph interactions are, sadly, lost when I export the Roassal visualizations to an image or even the javascript representation but as you can see in the screenshot from my Pharo image the number of infected on today's date is **297,236** individuals. Why is it so high? It seems as though our model is off by about 267,000 people!

We need to talk about testing. To date, the state of Florida has only tested 333,099 of it's 21.48 million inhabitants which is only 1.5 percent of it's population and, terrifyingly, of the 333,099 persons tested just over 9 percent tested positive for COVID-19. If we were to extrapolate that out to the total population then 1,988,663 people would currently test postive for the Corona virus. Extrapolating the ratio of infected to tested to the total population isn't going to give us an accurate estimation for a multitude of reasons but most notably because the 9% testing positive were presumably showing symptoms of the virus or are health care workers recently exposed to it. Nine percent of a total population seems quite high but [the CDC estimates that between 3-11% of the United States population are infected with the common flu each year](https://www.cdc.gov/flu/about/keyfacts.htm) and with the revalations that [carriers of the virus may remain asymptomatic](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7128959/) the **297,236** infected to date that our model shows might not be as far off as it seems.

Let's run our model on another population, one in which testing is more widely available and reporting presumably more accurate, New York City. First we do a little refactoring, DRYing up the code a bit by referencing the `populationSize` variable from the `state` array and changing the population size to match that of New York City according to [Google](https://www.google.com/search?rlz=1C1CHBF_enUS880US880&sxsrf=ALeKk01k46J52vuE8Es82o5UG-vbUeWQAA%3A1588005306283&ei=ugmnXrTlEKWMgge4g7O4Bg&q=population+of+new+york+city&oq=population+of+new+york+city&gs_lcp=CgZwc3ktYWIQAzIECCMQJzICCAAyAggAMgIIADICCAAyAggAMgYIABAHEB4yAggAMgIIADICCAA6BAgAEEdQ9KpFWNW2RWCTuEVoAHACeAGAAbYCiAHxEJIBCDAuMTEuMS4xmAEAoAEBqgEHZ3dzLXdpeg&sclient=psy-ab&ved=0ahUKEwi0qrrhhInpAhUlhuAKHbjBDGcQ4dUDCAw&uact=5). The [first reported case of Corona virus in New York City](https://www.wsj.com/articles/first-case-of-coronavirus-confirmed-in-new-york-state-11583111692) was also March 1, 2020 so our X-axis remains unchanged.

```smalltalk
| dt durationOfRecovery gamma populationSize r0 beta solver state system values stepper susceptible infected recovered b ds1 ds2 ds3 |

dt := 1.0.
durationOfRecovery := 7.5.
gamma := 1 / durationOfRecovery.
populationSize := 8399000.
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
state := { populationSize . 1 . 0 }.
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
ds1 interaction toggleDataset.
ds1 interaction popup text: [ :v | ((v key y) asInteger) asString, ' susceptible on ', ('March 1, 2020' asDate  + v key x day) asDate asString].
b add: ds1.

ds2 := RTData new.
ds2 label: 'Infected'.
ds2 noDot.
ds2 points: infected.
ds2 connectColor: Color red.
ds2 y: [ :v | v y ].
ds2 x: [ :v | v x ].
ds2 interaction toggleDataset.
ds2 interaction popup text: [ :v | ((v key y) asInteger) asString, ' infected on ', ('March 1, 2020' asDate  + v key x day) asDate asString].
b add: ds2.

ds3 := RTData new.
ds3 label: 'Recovered'.
ds3 noDot.
ds3 points: recovered.
ds3 connectColor: Color green.
ds3 y: [ :v | v y ].
ds3 x: [ :v | v x ].
ds3 interaction toggleDataset.
ds3 interaction popup text: [ :v | ((v key y) asInteger) asString, ' recovered on  ', ('March 1, 2020' asDate  + v key x day) asDate asString].
b add: ds3.

b addDecorator: (RTCursorFollower new color: Color gray).
b axisX
        title: '';
        labelRotation: -30;
        labelConversion: [ :v | ('March 1, 2020' asDate  + v day) asDate ].
b axisY title: 'Population'; noDecimal.
b legend right.

b build view canvas buildMorph extent: 1000@500; exportAsPNG.
```

![NYC COVID-19 transmission](https://i.imgur.com/YiHVSiH.png)

To review the accuracy of our NYC model we first reference the number of [Corona virus cases in NYC](https://www1.nyc.gov/site/doh/covid/covid-19-data.page) as of *April, 26, 2020* which is **153,204**. When we remove the Suscepible and Recovered lines and check the Infected on April, 26 we see our NYC model is much closer to reality with **147,108** infected individuals.

![NYC COVID-19 number of people infected](https://i.imgur.com/fzq2Wh2.png)

Now seems like a good time to abstract our code into classes to make it easier for us to run models. A live screencast is a great way to show how we can and turn our scripts into a collection of objects working together to achieve our goals with the the added bonus of extending the package to fit more models. I'll consider recording a live screencast within the near future.

In the meantime you can find the library I've already abstracted on Github at: [https://github.com/graves/2019-nCov](https://github.com/graves/2019-nCov)
