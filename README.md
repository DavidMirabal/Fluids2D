# Lax-Friedrichs 2D

Hecho por David Mirabal Betancort


## Uso del código:

Para ejecutar el código correctamente es necesario tener en el mismo directorio los siguientes archivos: `ci.py`, `dmb_figures.py`, `image_to_video.py`, `main.py`, `plot_mul_mom_angular.py`, `scheme.py` y `utils.py`

### Utilidad de cada acrhivo:

* `utils.py` : Clases y funciones útiles para cálculos. La clase *der* incluye el cálculo de la derivada euleriana en 2D como atributo. Se puede pedir la derivada con respecto a cualquiera de las dimensiones. También incluye el generador de trazadores (partículas que se posicionan sobre el fluido). Genera posiciones de N partículas a un radio dado. Contiene otra función que servirá para asignar el valor de velocidad del mapa de velocidades más cercano al trazador.

* `ci.py`: Contiene clases de condiciones iniciales. Cada clase es un tipo de condición inicial (ondas de sonido y vórtice Gresho). Llamando a funciones de la clase se pueden calcular las densidades, velocidades, ... iniciales.

* `dmb_figures.py`: Incluye una clase para homogeneizar las figuras de manera más sencilla y rápida. También incluye una función para generar colores a partir de un mapa de matplotlib.

> [!WARNING]
> Si hubiera algún problema para generar figuras con el código es probable que sea por el punto anterior. Hay que tener en cuenta que se usa LaTeX para el texto de las figuras y si no se tiene instalado LaTeX y los paquetes correspondientes no funcionará. En este caso se recomienda comentar todos los *plt.rcParams* de `dmb_figures.py`.

* `scheme.py`: Incluye el cálculo del método numérico Lax-Friedrichs para un tiempo dado. Llamando a las funciones de la clase se puede establecer las condicones de contorno periódicas y obtener las variables calculadas.

* `main.py`: Con este script se lanzará el código. Se lanza por terminal con un argumento 

    ```bash 
    $ python main.py config_file
    ```

    `config_file` es el archivo que definirá las propiedades de la simulación (veáse más adelante). Cuando se ejecuta generará frames (según se haya indicado en el archivo de configuración) con las densidades, velocidades, presión, momentos angulares, ... También generá archivos con los momentos angulares en función del tiempo.

* `image_to_video.py`: (Script creado para otras asignaturas) A partir de un directorio con los frames generados anteriormente, genera un video con las características deseadas. Se puede ejecutar tanto llamandolo por terminal como llamando a la clase que contiene. Por terminal se ejecuta así:     

    ```bash
    $ python image_to_video.py directorio_con_frames  fps directorio_salida titulo resize
    ```

    Se recomienda `resize=1` para mantener la calidad de los frames. Si los archivos de video empiezan a ser muy pesados si se recomienda bajar este valor. 

> [!WARNING]
>  Para ejecutar este último script es necesario tener instalado *FFmpeg*.

* `plot_mul_mom_angular.py`: Este script se ejecuta sin ningún argumento. Genera una figura donde muestra el momento angular en función del tiempo para multiples simulaciones. La variable *dir_stats* ha de definirse como el directorio donde se incluyen los csv que genera main.py de las simulaciones las cuales se quiere ver la evolución del momento angular. Además, intenta generar una figura de convergencia de la conservación del momento angular (para ello en el archivo de configuraciones se ha de escoger un nombre adecuado para la simulación, esto es poner NXXX donde XXX es el número de puntos por dimensión de la simulación).


### Archivo de configuración

Cada simulación se define por su archivo de configuración. En él se podrán cambiar parámetros de la visualización, de la simulación, de las condiciones iniciales y de los trazadores. Se dejan varios ejemplos de archivos de configuración que sirven como referencia. 

El archivo se divide en 5 bloques:

 * `DIRECTORIOS`: Incluye el nombre característico que se le quiera poner a la simulación, los directorios donde guardar los frames y donde guardar las evoluciones del momento angular, etc.

 * `VISUAL`: Incluye el número de frames que se quieres, el intervalo entre frames consecutivos y límites de los mapas de colores de los gráficos generados. También se puede modificar cada cuantos puntos se muestra una flecha representativa del vector velocidad en ese punto.

 * `SIMULACION`: Incluye el número de puntos en cada eje, los límites de la simulación, la densidad de fondo constante y el índice adiabático.

 * `CI`: Aquí se definirá si se quiere una simulación de ondas de sonido o si se quiere el vórtice. Se incluyen las características de las condiciones iniciales en el caso de que se usen las ondas de sonido.

 * `TRAZADORES`: Se definen las partículas trazadoras (su número y cuantas capas). Se generaran un número dado de partículas en radios que se escojan.


 ### Requerimientos
 
El código se ha testeado con Python 3.10, son obligatorios los siguientes paquetes externos a Python:

 * `Numpy`
 * `Matplotlib`

Como ya se ha mencionado, si se quiere generar video y no solo los frames se ha de instalar el software `FFmpeg`. Se usan los siguientes paquetes internos de Python: `os`, `sys`, `tqdm`, `configparser`, `subprocess` y `PIL`.

