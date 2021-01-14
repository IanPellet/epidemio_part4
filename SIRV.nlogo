turtles-own
  [ sick?                ;; if true, the turtle is infectious
    vaccinated?
    susceptible?
    remaining-recovered  ;; how many days of immunity the turtle has left
    vaccinated-time      ;; how long, in days, the turtle has been vaccinated
    sick-time            ;; how long, in days, the turtle has been infectious
    weakness
    age ]                ;; how many days old the turtle is

globals
  [ %susceptible         ;; what % of the population is susceptible
    %infected            ;; what % of the population is infectious
    %recovered           ;; what % of the population is recovered
    %vaccinated          ;; what % of the population is vaccinated
    lifespan             ;; the lifespan of a turtle
    mean-age
    chance-reproduce     ;; the probability of a turtle generating an offspring each tick
    carrying-capacity    ;; the number of turtles that can be in the world at one time
    recovered-duration   ;; how many days immunity lasts
    vaccine-duration     ;; how many days vaccine lasts
]

;; The setup is divided into four procedures
to setup
  clear-all
  setup-constants
  setup-turtles
  update-global-variables
  update-display
  reset-ticks
end

;; We create a variable number of turtles of which 10 are infectious,
;; and distribute them randomly
to setup-turtles
  create-turtles number-people
    [ setxy random-xcor random-ycor
      set age random-normal (lifespan / 2) sqrt(lifespan / 4)
      set sick-time 0
      set remaining-recovered 0
      set vaccinated-time 0
      set size 1.5  ;; easier to see
      set weakness random-float 10
      get-healthy ]
  ask n-of sick-people turtles
    [ get-sick ]
end

to get-sick ;; turtle procedure
  set sick? true
  set susceptible? false
  set vaccinated? false
  set remaining-recovered 0
  set vaccinated-time 0
end

to get-healthy ;; turtle procedure
  set susceptible? true
  set sick? false
  set vaccinated? false
  set remaining-recovered 0
  set sick-time 0
  set vaccinated-time 0
end

to become-recovered ;; turtle procedure
  set susceptible? false
  set sick? false
  set sick-time 0
  set remaining-recovered recovered-duration
end

to get-vaccine
  if age < max-age-for-vaccine * 365 [
    if sick? = false [
      if random-float 1000 < ( (vaccine-chance / 10 ) * 100 ) [
         set susceptible? false
         set sick? false
         set vaccinated? true
         set sick-time 0
         set vaccinated-time 0
      ]
    ]
  ]
end

;; This sets up basic constants of the model.
to setup-constants
  set lifespan 80 * 365     ;; 80 times 365 days = 80 years
  set carrying-capacity 300
  set chance-reproduce 1
  set recovered-duration 5 * 365 ;; 5 times 365 days = 5 years
  set vaccine-duration vaccine-efficiency * 365  ;; 10 times 365 days = 10 years
end

to go
  if all? turtles [ sick? = false]
    [ stop ]

  if all? turtles [ sick? = true]
    [ stop ]

  ask turtles with [sick? = true]
    [ infect ]

  ask turtles with [sick? = false]
    [ if age > (15 * 365) and age < (60 * 365) [ reproduce ] ]

  ask turtles with [sick-time >= duration]
    [ifelse weakness > 2.5 [recover-or-die][get-sick]]

  ask turtles with [susceptible? = true]
   [ if vaccine? [get-vaccine] ]

  ask turtles with [recovered?]
    [ if remaining-recovered = 1 [ set susceptible? true] ]

  ask turtles with [vaccinated-time >= vaccine-duration]
    [get-healthy]

  ask turtles with [ age > (70 * 365)]
      [ maybe-die ]

  ask turtles [
   get-older
   move
  ]

  update-global-variables
  update-display
  tick
end


to maybe-die
  random-seed 23
  if random-float 1 < (age / lifespan)
    [ die ]
end

to update-global-variables
  if count turtles > 0
    [ set %susceptible (count turtles with [ susceptible? ] / count turtles) * 100
      set %infected (count turtles with [ sick? ] / count turtles) * 100
      set %recovered (count turtles with [ recovered? ] / count turtles) * 100
      set %vaccinated (count turtles with [ vaccinated? = true ] / count turtles) * 100
      set mean-age mean [age / 365] of turtles
  ]
end

to update-display
  ask turtles[
      if shape != turtle-shape [ set shape turtle-shape ]
      if vaccinated? = true [set color yellow ]
    set color ifelse-value sick? [ red ][ ifelse-value vaccinated? = true [ yellow ] [ ifelse-value recovered? [blue] [green] ] ]
  ]
end

;;Turtle counting variables are advanced.
to get-older ;; turtle procedure
  ;; Turtles die of old age once their age exceeds the
  ;; lifespan (set at 50 years in this model).
  set age age + 1
  if age > lifespan [ die ]
  if recovered? [ set remaining-recovered remaining-recovered - 1]
  if vaccinated? = true [ set vaccinated-time vaccinated-time + 1]
  if sick? [ set sick-time sick-time + 1 ]
end

;; Turtles move about at random.
to move ;; turtle procedure
  rt random 100
  lt random 100
  fd 1
end

;; If a turtle is sick, it infects other turtles on the same patch.
;; Immune and vaccinated turtles don't get sick.
to infect ;; turtle procedure
  ask other turtles-here with [ susceptible? = true]
    [ if random-float 100 < infectiousness
      [ get-sick ] ]
end


;; Once the turtle has been sick long enough, it
;; either recovers (and becomes immune) or it dies.
to recover-or-die                                ;; turtle procedure
  if sick-time > duration                        ;; If the turtle has survived past the virus' duration, then
    [ ifelse random-float 100 < chance-recover   ;; either recover or die
      [ become-recovered ]
      [ die ] ]
end

;; If there are less turtles than the carrying-capacity
;; then turtles can reproduce.
to reproduce
  if count turtles < carrying-capacity and random-float 5000 < chance-reproduce
    [ hatch 1
      [ set age 1
        lt 45 fd 1
        get-healthy ] ]
end


to-report recovered?
  report remaining-recovered > 0
end


to startup
  setup-constants ;; so that carrying-capacity can be used as upper bound of number-people slider
end
@#$#@#$#@
GRAPHICS-WINDOW
455
86
1158
640
-1
-1
10.7
1
10
1
1
1
0
1
1
1
-32
32
-25
25
1
1
1
ticks
45.0

SLIDER
15
95
220
128
duration
duration
0.0
99.0
14.0
1.0
1
days
HORIZONTAL

SLIDER
235
10
429
43
chance-recover
chance-recover
0.0
99.0
85.0
1.0
1
%
HORIZONTAL

SLIDER
235
50
430
83
infectiousness
infectiousness
0.0
99.0
70.0
1.0
1
%
HORIZONTAL

BUTTON
15
245
85
280
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
91
245
162
281
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

PLOT
15
340
430
640
Populations
days
people
0.0
365.0
0.0
200.0
true
true
"" ""
PENS
"Susceptible" 1.0 0 -13840069 true "" "plot count turtles with [ not sick? and not recovered?]"
"Infectious" 1.0 0 -2674135 true "" "plot count turtles with [ sick? ]"
"Recovered" 1.0 0 -13791810 true "" "plot count turtles with [ recovered? ]"
"Vaccinated" 1.0 0 -1184463 true "" "plot count turtles with [ vaccinated? = true ]"

SLIDER
15
8
220
41
number-people
number-people
10
carrying-capacity
120.0
1
1
NIL
HORIZONTAL

MONITOR
550
25
635
70
NIL
%infected
3
1
11

MONITOR
835
25
925
70
years
ticks / 365
3
1
11

MONITOR
635
25
720
70
NIL
%recovered
3
1
11

MONITOR
720
25
812
70
%vaccinated
%vaccinated
3
1
11

CHOOSER
235
190
430
235
turtle-shape
turtle-shape
"person" "circle" "turtle"
0

SLIDER
15
50
220
83
sick-people
sick-people
0
number-people
10.0
1
1
NIL
HORIZONTAL

SLIDER
235
95
430
128
vaccine-efficiency
vaccine-efficiency
0
25
6.0
1
1
years
HORIZONTAL

SLIDER
235
140
430
173
vaccine-chance
vaccine-chance
0
100
20.0
1
1
%
HORIZONTAL

SLIDER
15
140
220
173
max-age-for-vaccine
max-age-for-vaccine
0
lifespan / 365
10.0
1
1
years
HORIZONTAL

SWITCH
15
200
125
233
vaccine?
vaccine?
0
1
-1000

MONITOR
455
25
550
70
NIL
%susceptible
3
1
11

MONITOR
925
25
1002
70
NIL
mean-age
3
1
11

MONITOR
1000
25
1072
70
Total pop
count turtles
3
1
11

@#$#@#$#@
\begin{itemize}
    \item [$\bullet$] \textbf{remaining-recovered} : Représente le nombre de jours passés depuis la guérison de la personne. Pendant cette période où la personne est guérie du virus, elle a une certaine immunité contre celui-ci. Cette immunité n'est pas définitive et dès que cette période atteint son maximum ( variable \textbf{recovered-duration} ) la personne redevient susceptible ;
    \item [$\bullet$] \textbf{sick-time} : Représente la durée depuis laquelle la personne a contracté le virus. Cette durée s'exprime en jour, commence à $0$ et finie lorsque \textbf{sick-time $=$ duration} ;
    \item [$\bullet$] \textbf{duration} : Représente la durée pendant laquelle une personne est malade et infectieuse, sa valeurs est arbitraire car laissée au choix de l'utilisateur via le curseur;
    \item [$\bullet$] \textbf{vaccinated-time} : Représente le nombre de jours passés depuis que la personne à reçu le vaccin. Pendant cette période où la personne invulnérable du virus car vaccinée. Cette invulnérabilité n'est pas définitive et dès que cette période atteint son maximum ( variable \textbf{vaccine-duration } ) la personne redevient susceptible ;
    \item [$\bullet$] \textbf{vaccinae-duration} : Cette variable représente la durée en années pendant laquelle le vaccin est efficace contre le virus, sa valeur est laissée au choix à l'utilisateur dans l'interface via le curseur \textbf{vaccine-efficiency} ; 
    \item [$\bullet$] \textbf{carrying-capacity} : le nombre d'individus qui peuvent se trouver dans le monde à un moment donné, lorsqu'elle est atteinte les individus ne peuvent plus se reproduire. ;
    \item [$\bullet$] \textbf{lifespan } : Représente l'espérance de vie des individus, elle est fixée à $80$ ans ;
    \item [$\bullet$] \textbf{chance-reproduce } : Probabilité qu'une tortue se reproduise à chaque tour, fixée à $\frac{1}{1000}$ (au delà on observe une reproduction dégénérée des individus) ;
    \item [$\bullet$] \textbf{weakness} : Désigne la fragilité d'un individu. C'est un flottant (entre 1 et 10), définie aléatoirement pour chaque individu. Il accentue le fait que certains individus sont plus vulnérables que d'autres et que le temps qu'ils doivent mettre pour guérir du virus sera plus long que la normale ( définie sur \textbf{duration}); 
    \item [$\bullet$] \textbf{max-age-for-vaccine} : Désigne l'âge maximal après lequel il n'est plus possible de recevoir le vaccin. Sa valeur en années est laissée au choix à l'utilisateur; 
    \item [$\bullet$] \textbf{vaccine-chance} : Cette variable est peut-être la variable la plus "bricolée" du modèle, elle permet de maintenir une certaine stabilité dans le comportement de celui-ci. Biologiquement elle représenterait le taux d'efficacité du vaccin ainsi que le pourcentage de la population ayant reçu le vaccin. Elle est laissé au choix à l'utilisateur, mais par défaut il vaut mieux la laisser sur une valeur faible ( comme $2\%$) pour ne pas dégénérer le modèle. 
    \item [$\bullet$] \textbf{infectiousness} : Probabilité d'infecter un autre individu en cas de contact
    \item [$\bullet$] \textbf{change-recover} : Une fois la période infectieuse passée il y $3$ possibilité : 
        \begin{enumerate}
            \item Si la variable \textbf{weakness} de la personne est trop faible, celle-ci redevient malade à nouveau ( sa période infectieuse est multipliée par $2$ ;
            \item Si la probabilité \textbf{chance-recover} est assez élevée la personne devient saine et guérie avec une période d'immunité \textbf{recovered-duration} ;
            \item Sinon elle meurt. 
        \end{enumerate}
\end{itemize}s
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.1.1
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
1
@#$#@#$#@
