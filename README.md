# single-chamber-decompression
This code is an implementation of the solution for the single chamber rapid decompression problem, described in the article: "Gasdynamics of rapid and explosive decompressions of pressurized aircraft including active venting" by Alfonso Pagani and Erasmo Carrera

OBS: IMPORTANT! The Equation (2) and (5) are wrongly written in the article. The correct ones:

Pi*/Pj = ((gamma + 1)/2)^(gamma/(gamma - 1))

mij_dot = Aeff * sqrt(2 * Pi * rhoi * (gamma/(gamma - 1))*((Pj/Pi)(2/gamma)) * (1 - (Pj/Pi)((gamma - 1)/gamma)))

Although this mistake in the article, the results are correct.
