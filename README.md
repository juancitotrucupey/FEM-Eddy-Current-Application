# FEM-Eddy-Current-Application
Se aplica el método de elemento finitos para resolver las ecuaciones magnéticas cuasi-estáticas (eddy current approximation of Maxwell equations)

El material de este trabajo fue tomado completamente del repositorio https://github.com/rosskynch/MIT_Forward creado por Ross Kynch.
Se tomarón los códigos creados por Ross Kynch y se estudió la implementación del problema benchmark 6 del TEAM:
Un cascarón (en este caso una esfera completa) esférico de conductor dentro de un campo magnético uniforme que 
oscila armónicamente. Se calcula la corriente de Foucault dentro del conductor y se calcula el campo electromagnético.

El archivo informe.pdf contiene una breve exposición del problema, los resultados de un par de simulaciones y una explicación
superficial de la implementación del código, el cual fue creado a partir del deal.ii.

Si se desea compilar y ejecutar el código en este repositorio se recomeinda instalar deal.ii versión 8.3 actualiza al
commit (SHA hash 79583e56..) del 6 de julio

El archivo instruccione.md contiene los comandos que se deben utilizar para compilar y ejecutar, al igual que una breve exposición
(realizada por Ross Kynch) sobre el problema original que motivo la creación de este software.

Se recomienda consultar el siguiente artículo, el cual expone la forma en la cual se implementaron los elementos de Nedeléc en
hexaedros y se utilizan dentro del código:

R.M. Kynch, P.D. Ledger, "Resolving the sign conflict problem for hp–hexahedral Nédélec elements
with application to eddy current problems",Computers and Structures 181 (2017) 41–54, http://dx.doi.org/10.1016/j.compstruc.2016.05.021

