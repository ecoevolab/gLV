#
#
# library(markovchain)
#
# transition_matrix <- t( apply(output, 1, function(row) row / sum(row)) )
#
# specsNames = paste0("State", 1:nrow(output))
# statesNames <- paste0("State", 1:ncol(output))
#
# markovB <- new("markovchain", states = statesNames,
#                transitionMatrix = transition_matrix,
#                dimnames = list(specsNames, statesNames),
#                name = "A markovchain Object"
# )
#
# steadyStates(markovB)
#
# #-----------------------------------------------------------------------------------------------------------#
# library(rootSolve)
#
# # Definir la función del sistema de ecuaciones
# lv_equations <- function(t, state, params) {
#
#   with (as.list(c(state,params)), {
#   r <- params$Growths         # Tasas de crecimiento
#   A <- params$Interactions    # Matriz de interacciones
#
#   dxdt <- x * (r + A %*% x)   # dx/dt = x * (r + A * x)
#   return(list(dxdt))
#   })
# }
#
# # Condiciones iniciales
# x_init <- params$Population  # Poblaciones iniciales
# # x_init <- c(.5, .5)
#
# # Resolver los estados estacionarios
# steady_state <- steady(y = x_init,
#                        func = lv_equations,
#                        parms = params,
#                        method = "runsteady")
#
# # Mostrar los estados estacionarios
# print(steady_state$y)
#
#
#
# #-----------------------------------------------------------------------------------------------------------#
# # Encontrar la Jacobiana
# library(numDeriv)
#
# r <- params$Growths         # Tasas de crecimiento
# A <- params$Interactions    # Matriz de interacciones
#
# # Define la función que representa el sistema de ecuaciones
# lv_equations <- function(x) {
#   n <- length(x)
#   f <- numeric(n)
#   for (i in 1:n) {
#     f[i] <- x[i] * (r[i] + sum(A[i, ] * x))
#   }
#   return(f)
# }
#
# # Define el punto en el que quieres evaluar la Jacobiana
# x <- params$Population  # Poblaciones iniciales
#
# # Calcula la matriz Jacobiana
# J <- jacobian(lv_equations, x)
# print(J)
#
# # Calculate eigenvalues
# eigenvalues <- eigen(J)$values
# print(eigenvalues)
#
# # Analyze stability based on eigenvalues
# if (all(Re(eigenvalues) < 0)) { # ALL are true
#   print("The steady state is locally stable.")
# } else if (any(Re(eigenvalues) > 0)) { # One is not TRUE
#   print("The steady state is unstable.")
# } else {
#   print("The steady state is a saddle point.")
# }


