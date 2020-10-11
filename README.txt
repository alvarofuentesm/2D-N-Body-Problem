Se proporcionan unos ejemplos de muestra en los archivos
GIF-GPU-Resultado-20000x20 iteraciones.gif
GIF-GPU-Resultado-10000 iteraciones.gif
GIF-CPU-Resultado-10000 iteraciones.gif

Intrucciones de uso:

- Asegurar la existencia de las carpetas /data y /gif
- Compilar con nvcc n_body.cu
- Ejecutar .exe generado

Se ejecutaran automáticamente los tres métodos usando como input el archivo input.txt de 500 cuerpos.

- Para visualizar resultados abrir el notebook Visualizador_N_Body.ipynb 
y ejecutar celdas bajo los títulos de: 

Importar librerias
Plot and generate GIF functions
Examples

NOTA IMPORTANTE: En caso de tener problemas con alguna librería en el visualizador lo mejor es abrir el notebook en 
Google Colab y subir las carpeta /data y /gif

- Los GIF generados se encuentran en la carpeta /gif

Nota: No se genera un GIF para la implementación basada de Fast N-Body.
Nota 2: Es posible generar un nuevo archivo input.txt con la celda Random Generator
Nota 3: Los cuerpos representados en las imágenes/GIF son proporcionales a su masa, no su diamétro, pero en la simulación son 
tratados como partículas puntuales que se pueden atravesar.  