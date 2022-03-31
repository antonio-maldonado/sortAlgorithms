// compilacion y ejecuta :     g++ -Wall testsort0.cpp -o sort -lrt && ./sort

// define encabezdos de funciones a utlizar
#include<iostream>
#include<cstdio>
#include<cstdlib>
#include<ctime>
#include<cmath>
#include <unistd.h>
using namespace std;

// *************************************************************
//    funciones para despliegue de datos
// *************************************************************
void muestraUso( string );
void DespliegaValores(FILE*,int*,int,int,const char*);
FILE *gnuplot_pipe;
int totalDatos = 0;
int counter = 0;

// *************************************************************
//    definicion de funciones comunes y algoritmos de ordenamiento
// *************************************************************
typedef int Item;
#define key( A ) ( A )
#define less( A, B ) ( key(A) < key(B) )
#define exch( A,B ) { Item t = A; A = B; B = t; } 
#define compexch( A, B ) if ( less(B,A) ) exch(A,B)
#define min( A, B ) (A < B) ? A : B
#define eq( A, B ) ( !less(A,B) && !less(B,A) )
#define M 10

// *************************************************************
//    Funcion seleccion metodo
// *************************************************************

//Metodo Selection 
void selection(Item a[], int l, int r)
  { int i, j;
    for (i = l; i < r; i++)
      { int min = i;
        for (j = i+1; j <= r; j++) 
            if (less(a[j], a[min])) min = j;
        exch(a[i], a[min]);
         DespliegaValores( gnuplot_pipe, a, totalDatos, ++counter, "Selection sort" );
      } 
  }

//Metodo Insertion
void insertion(Item a[], int l, int r)
  { int i;
    for (i = l+1; i <= r; i++) compexch(a[l], a[i]);
    for (i = l+2; i <= r; i++)
      { int j = i; Item v = a[i]; 
        while (less(v, a[j-1]))
          { a[j] = a[j-1]; j--; }
        a[j] = v; 
DespliegaValores( gnuplot_pipe, a, totalDatos, ++counter, "Insertion Sort" );
     } 
  }
//Metodo Insetion para el quick
void insertion(Item a[], int l, int r, bool flag)
  { int i;
    for (i = l+1; i <= r; i++) compexch(a[l], a[i]);
    for (i = l+2; i <= r; i++)
      { int j = i; Item v = a[i]; 
        while (less(v, a[j-1]))
          { a[j] = a[j-1]; j--; }
        a[j] = v; 
        if(flag)
        DespliegaValores( gnuplot_pipe, a, totalDatos, ++counter, "insertion sort" );
      } 
  }
//Metodo Metodo Bubble  
void bubble(Item a[], int l, int r)
  { int i, j;
    for (i = l; i < r; i++){
      for (j = r; j > i; j--)
        compexch(a[j-1], a[j]);
         DespliegaValores( gnuplot_pipe, a, totalDatos, ++counter, "bubble sort" );
         }
  } 
  
//Metodo shell knuth https://en.wikipedia.org/wiki/Shellsort
 
void shellsort(Item a[], int l, int r) 
  { 
  int i, h; 
    for (h = 1; h <= (r-l)/9; h = 3*h+1) ;
    for ( ;h > 0; h /= 3)
      for (i = l+h; i <= r; i++)
        { int j = i; Item v = a[i]; 
          while (j >= l+h && less(v, a[j-h]))
            { a[j] = a[j-h]; j -= h; }
          a[j] = v; 
          DespliegaValores( gnuplot_pipe, a, totalDatos, ++counter, "Shell Knuth" );
        } 
  }
  
//Metodo Shell 1959
void shell(Item a[], int l, int r) 
  { 
  int i, h, k; 
      for (k=1, h=(r-l+1)/2; h>0; h=(r-l+1)/int(pow(2,++k)))
      for (i = l+h; i <= r; i++)
        { int j = i; Item v = a[i]; 
          while (j >= l+h && less(v, a[j-h]))
            { a[j] = a[j-h]; j -= h; }
          a[j] = v; 
          DespliegaValores( gnuplot_pipe, a, totalDatos, ++counter, "Shell 1959" );
        } 
  }
 
//Metodo Frank & Lazarus, 1960
void shellF(Item a[], int l, int r) 
  { 
  int i, h, k; 
      for (k=1, h=1+(2*(r-l+1)/4); h>0; h=2*(r-l+1)/int(pow(2,(++k)+1)))
      for (i = l+h; i <= r; i++)
        { int j = i; Item v = a[i]; 
          while (j >= l+h && less(v, a[j-h]))
            { a[j] = a[j-h]; j -= h; }
          a[j] = v; 
          DespliegaValores( gnuplot_pipe, a, totalDatos, ++counter, "Shell Frank y Lazarus, 1960" );
        } 
  }
  
  int partition(Item a[], int l, int r);

//Metodo Quick Básico
void quicksort(Item a[], int l, int r)
  { int i;
    if (r <= l) return;
    i = partition(a, l, r);
    quicksort(a, l, i-1);
    quicksort(a, i+1, r);
  }
  int partition(Item a[], int l, int r)
  { int i = l-1, j = r; Item v = a[r];
    for (;;)
      { 
        while (less(a[++i], v)) ;
        while (less(v, a[--j])) if (j == l) break;
        if (i >= j) break;
        exch(a[i], a[j]);
      }
    exch(a[i], a[r]);
    DespliegaValores( gnuplot_pipe, a, totalDatos, ++counter, "Quick Básico" );
    return i;
  }
//Metodo Quick Insertion
  void quicksortM(Item a[], int l, int r)
  { int i; 
    if (r-l <= M) return;
    exch(a[(l+r)/2], a[r-1]);
   compexch(a[l], a[r-1]); 
      compexch(a[l], a[r]); 
        compexch(a[r-1], a[r]);
    i = partition(a, l+1, r-1);
    quicksortM(a, l, i-1);
    quicksortM(a, i+1, r);
  }
  
//Metodo Quick Mediana
void quicksort3(Item a[], int l, int r)
  { int i, j, k, p, q; Item v;
    if (r <= l) return;
    v = a[r]; i = l-1; j = r; p = l-1; q = r;
    for (;;)
      { 
        while (less(a[++i], v)) ;
        while (less(v, a[--j])) if (j == l) break;
        if (i >= j) break;
        exch(a[i], a[j]);
        if (eq(a[i], v)) { p++; exch(a[p], a[i]); }
        if (eq(v, a[j])) { q--; exch(a[q], a[j]); }
      }
    exch(a[i], a[r]); j = i-1; i = i+1;
    for (k = l  ; k < p; k++, j--) exch(a[k], a[j]);
    for (k = r-1; k > q; k--, i++) exch(a[k], a[i]);
    DespliegaValores( gnuplot_pipe, a, totalDatos, ++counter, "Quick Mediana" );
    quicksort3(a, l, j);
    quicksort3(a, i, r); 
  }

//Metodo insertion para Merge
void insertion2(Item a[], int l, int r)
{ 
  for (int i = l+1; i <= r; i++) 
    compexch(a[l], a[i]);
  for (int i = l+2; i <= r; i++)
   { 
     int j = i; Item v = a[i]; 
       while (less(v, a[j-1]))
        { a[j] = a[j-1]; j--; }
         a[j] = v;
     DespliegaValores( gnuplot_pipe, a, totalDatos, ++counter, "Insertion2 Sort" );
     } 
  };

//Metodo Merge Sort Para Uso Despues  

Item *aux = NULL;
void merge(Item a[], int l, int m, int r)
  { int i, j, k;

    for (i = m+1; i > l; i--) 
    aux[i-1] = a[i-1];
    for (j = m; j < r; j++) 
    aux[r+m-j] = a[j+1];
    for (k = l; k <= r; k++)
       if (less(aux[i], aux[j])) 
          a[k] = aux[i++]; 
        else a[k] = aux[j--];
DespliegaValores( gnuplot_pipe, a, totalDatos,++counter, "Merge Sort" );
  }


//Metodo Merge Basico Top Down 

void mergesort(Item a[], int l, int r)    
  { int m = (r+l)/2;
    if (r <= l) 
    return;
    mergesort(a, l, m);  
    mergesort(a, m+1, r);
    merge(a, l, m, r);

  }


//Metodo Merge Basico Mejorado  

void mergeAB(Item c[], Item a[], int N, Item b[], int L )
  { int i, j, k;
    for (i = 0, j = 0, k = 0; k < N+L; k++)
      {
        if (i == N) { c[k] = b[j++]; continue; }
        if (j == L) { c[k] = a[i++]; continue; }
        c[k] = (less(a[i], b[j])) ? a[i++] : b[j++];
      }
DespliegaValores( gnuplot_pipe, a, totalDatos,++counter, "MergeAB Sort" );
  }


void mergesortABr(Item a[], Item b[], int l, int r)
  { int m = (l+r)/2;
    if (r-l <= 10) { insertion(a, l, r); return; }
    mergesortABr(b, a, l, m);  
    mergesortABr(b, a, m+1, r);
    mergeAB(a+l, b+l, m-l+1, b+m+1, r-m);
  }


void mergesortAB(Item a[], int l, int r)
  { int i;
    for (i = l; i <= r; i++) 
    aux[i] = a[i];
    mergesortABr(a, aux, l, r);
    DespliegaValores( gnuplot_pipe, a, totalDatos,++counter, "Merge Sort" );  
}



//Metodo Merge Sort Bottom-up 
#define min(A, B) (A < B) ? A : B
void mergesortBU(Item a[], int l, int r)
  { int i, m;
    for (m = 1; m < r-l; m = m+m)
      for (i = l; i <= r-m; i += m+m)
        merge(a, i, i+m-1, min(i+m+m-1, r));

  }

#define pq(A) a[l-1+A]
void fixDown(Item a[], int k, int N)
  { int j;
    while (2*k <= N)
      { j = 2*k;
        if (j < N && less(a[j], a[j+1])) j++;
        if (!less(a[k], a[j])) break;
        exch(a[k], a[j]); k = j;
        DespliegaValores( gnuplot_pipe, a, totalDatos,++counter, "Heap Sort" ); 
     }

  }
void heapsort(Item a[], int l, int r)
  
{ 
int k, N= r-l+1;
    for (k = N/2; k >= 1; k--) 
      fixDown(&pq(0), k, N);
    while (N > 1) 
      { exch(pq(1), pq(N)); 
        fixDown(&pq(0), 1, --N); }
  }

void sort( Item *a, int l, int r, int S )
{
  // inicializa variables
  totalDatos = r - l + 1;
  DespliegaValores( gnuplot_pipe, a, totalDatos, counter, "Datos Iniciales" ); 

  // selecciona algoritmo de ordenamiento
  switch( S )
    {
      case 1:
        cout << "\tSelection sort" << endl;
        selection(a, l, r);
        break;
        
        case 2:
        cout << "\tInsertion sort" << endl; //Comprartido con shell frank
        insertion(a, l, r,true);
        break;
        
        case 3:
        cout << "\tBubble sort" << endl;
        bubble(a, l, r);
        break;
        
        case 4:
        cout << "\tShell Knuth" << endl;
        shellsort(a, l, r);
        break;
        
        case 5:
        cout << "\tShell 1959" << endl;
        shell(a, l, r);
        break;

        case 6:
        cout << "\tQuick Básico" << endl;
        quicksort(a, l, r);
        break;
        
       case 7:
       cout << "\tQuick Insertion" << endl;
       quicksortM(a, l, r),  insertion(a, l, r,false); 
        break;
       case 8:
       cout << "\tQuick  Mediana" << endl;
       quicksort3(a, l, r);
        break;
       
        case 9: 
	aux = new Item[totalDatos];
	mergesort(a,l,r);
	cout<< "\n" << "Metodo Merge Sort Basico"<<endl;
        delete []aux;
	break;

        case 10: 
	aux=new Item[totalDatos];
	mergesortAB(a,l,r);
	cout<< "\n" << "Metodo Merge Basico Mejorado Con Insertion"<<endl;
        delete[]aux;
	break;

        case 11: 
	aux=new Item[totalDatos];
	mergesortBU(a,l,r);
	cout<< "\n" << "Metodo Merge Sort Bottom Up "<<endl;
        delete[]aux;
	break;
        
        case 12: 
	cout<< "\n" << "Metodo Heap Sort "<<endl;
        heapsort(a,l,r);
        break;
        
        default:
        cout << " No existe el metodo seleccionado..." << endl;
        break;  
    }
}        



// *************************************************************
//    Funcion principal
// *************************************************************
int main( int argc, char *argv[] )
{
  // declara variables a utilizar
  timespec ts_beg, ts_end;          // estructuras de tiempo
  int N = 1000;                             // numero de datos  
  int Selec = 1;                             // metodo de ordenamiento

  // valores de la terminal 
  if ( argc > 1 )
    {
      // selecciona opciones
      switch( argc )
        {
          case 3:
            N = (int) atoi( argv[1] );
            Selec = (int) atoi( argv[2] );
            break;
        }
    }
  else
    {
      muestraUso( argv[0] );
      return EXIT_SUCCESS;
    }
        

  
  // inicializa variables de memoria y visualizacion
  Item *a = new Item[N];             // separa memoria
  gnuplot_pipe = popen( "gnuplot", "w" );

  // genera valores aleatorios
  for ( int i = 0; i < N; i++ )
    a[i] = rand()%1000;

  // ordena valores
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts_beg);
  sort( a, 0, N-1, Selec );
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts_end);
                      
  // despliega resultados
  float tiempo = (float(ts_end.tv_sec) - float(ts_beg.tv_sec))*1.0e3
                      + (float(ts_end.tv_nsec) - float(ts_beg.tv_nsec))/1.0e6;
  cout << "Tiempo empleado (msec) : " << tiempo << endl;
  sleep( 1 );


  // termina programa
  pclose( gnuplot_pipe );  
  delete [] a;
  return EXIT_SUCCESS;
}



// *************************************************************
//    Funcion muestra forma de uso del programa 
// *************************************************************
void muestraUso( string name )
{
    std::cerr << "Uso: " << name << " Numero_Elementos Metodo" << endl
              << "\tMetodo: " << endl
              << "\t\t1. Selection Sort" << endl
              << "\t\t2. Insertion Sort" << endl
              << "\t\t3. Bubble Sort" << endl
              << "\t\t4. Shell Knuth" << endl
              << "\t\t5. Shell 1959" << endl
              << "\t\t6. Quick Básico" << endl
              << "\t\t7. Quick Insertion" << endl
              << "\t\t8. Quick  Mediana" << endl
              << "\t\t9. Metodo Merge Sort Basico" << endl
              << "\t\t10. Metodo Merge Basico Mejorado Con Insertion" << endl
              << "\t\t11. Metodo Merge Sort Bottom Up" << endl
              << "\t\t12. Heap Sort" << endl
//              << "\tOpciones: " << endl
              << endl              
              << "\tEjemplo: " << name << " 10000 1"
              << endl;
}



// *************************************************************
//    Funcion visulizacion
// *************************************************************
void DespliegaValores( FILE* salida, int* V, int numDatos, int step, const char* title ) //valores de las imagenes
{
  // opciones de trabajo
  bool flagG = false;    // true si se quiere graficos, para saber si quiero saber el tiempo de ordenamiento sin graficos poner false
  bool flagF = false;    // true si se quiere archivos(imagenes), false salida en ventana

  // ciclo de despliegue  Â¿BUG del GNUPLOT?
  for ( int k = 0; k < 2; k++ )
    {
      // salida de los graficos?
      if ( flagG )
        {
          // salida en modo grafico o en archivo
          if ( flagF )  // archivo !!!
          {
            fprintf( salida, "set terminal jpeg enhanced\n");
            fprintf( salida, "set output \"imagenes/out%04i.jpg\"\n", step );
          }
        fprintf( salida, "set title \"%s\"\n", title );
        fprintf( salida, "set xlabel \"orden\"\n" );
        fprintf( salida, "set ylabel \"valor\"\n" );
        fprintf( salida, "set xrange [%f:%f]\n", 0.0, float(numDatos) );
        fprintf( salida, "set yrange [%f:%f]\n", 0.0, 1000.0 );
        fprintf( salida, "plot '-' using 1:2 title '' with points pointtype 6\n" );
        for ( int c = 0; c < numDatos; c++ )
          fprintf( salida, "%f %f\n", float(c), float(V[c]) );
        fprintf( salida, "e\n" );
        fflush( salida );
        usleep( 10000);
      }
    }  
}
