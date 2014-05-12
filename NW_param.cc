// switch off assertions
// #define NDEBUG
#include <assert.h>
#include <ctype.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "NW_type.h"
#include "NW_misc.h"
#include "NW_param.h"


UNCHECKED_PARAM_TABLE::UNCHECKED_PARAM_TABLE(const char* _file_name)
{
  head = tail = NULL; 
  default_param_name = "first_params";
  current_param = strdup(default_param_name); 
  default_log = NULL;
  if (_file_name) {
    file_name = strdup(_file_name);
    parse_file(file_name);
  } else
    file_name = NULL;
}

UNCHECKED_PARAM_TABLE::~UNCHECKED_PARAM_TABLE() 
{ 
  if (default_log)
    fclose(default_log); 
  if (file_name)
    delete file_name;
}

PARAM *
UNCHECKED_PARAM_TABLE::put(const char* name)
{
  PARAM *pe;

  pe = find(name);
  if (pe) 
    return(pe);
  pe = new PARAM(strdup(name));;
  pe->size = 0;
  pe->data = NULL;
  if (head == NULL) {
    head = tail = pe;
    pe->next = NULL;
  } else {
    pe->next = head;
    head = pe;
  }
  return(pe);
}


PARAM *
UNCHECKED_PARAM_TABLE::find(const char* name)
{
  PARAM *pe;

  for(pe = head; pe != NULL; pe = pe->next) {
    if (strcmp(pe->name,name) == 0)
      return(pe);
  }
  return(NULL);
}

void
UNCHECKED_PARAM_TABLE::del(const char* )
{
// FIX
}

void
PARAM::resize(int new_size)
{
  char **ndata;

  if (new_size <= size)
    return;
  ndata = new CHAR*[new_size];
  if (data != NULL) {
    memcpy(ndata,data,size * sizeof(char*));
    delete[] data;
  }
  data = ndata;
  size = new_size;
}

void
UNCHECKED_PARAM_TABLE::put(const char* name, int index, const char* string) 
{
  PARAM *pe;
  
  pe = put(name);
  pe->resize(index+1);
  pe->data[index] = strdup(string);
}

void
UNCHECKED_PARAM_TABLE::parse_token(const char* token, PARAM_MODE mode)
{
// check for -switchname
  if (token[0] == '-'
      && !(isdigit(token[1]) || token[1] == '.')) {

    // setup new parameter name = switchname
    delete current_param;

    if (find(&token[1]) == NULL) {
      // parameter does not exist
      current_param = strdup(&token[1]);
      put(current_param);
    } else {
      // parameter already exists
      if (mode == UNDERRIDE)
	current_param = NULL;
      else {
	if (mode == OVERRIDE) {
	  del(&token[1]);
	}
	// else mode == APPEND
	current_param = strdup(&token[1]);
      }
    }

  } else {

    // put another value onto current parameter (if not ignoring it)
    if (current_param != NULL)
      put(current_param, get_size(current_param), token);
  }
}

void
UNCHECKED_PARAM_TABLE::parse_cmd_line(INT argc, const char* *argv, PARAM_MODE mode)
{
  INT i;
  start_parse();
  for(i=1; i<argc; i++) {
    parse_token(argv[i], mode);
  }
}

void
UNCHECKED_PARAM_TABLE::parse_file(const char* filename, PARAM_MODE mode)
{
  FILE* file;
  const int maxline = 2000;
  char* line = new CHAR[maxline];

  start_parse();
  file_name = strdup(filename);
  file = fopen((char*)filename,"r");
  if (file == NULL) {
    perror(filename);
    exit(-1);
  }
  while(fgets(line, maxline, file) != NULL) {
    if (parse_line(line, mode))
      break;
  }
  fclose(file);
  delete line;
}

BOOL
UNCHECKED_PARAM_TABLE::parse_line(CHAR* ptr, PARAM_MODE mode)
{
  CHAR* token;

  // check for non-ascii chars
  for(token=ptr; *token; token++) {
    if (!(isprint(*token) || isspace(*token))) {
      panic("bad char in pfile = 0x%x, in following line:\n%s\n\n", *token,
	    ptr);
    }
  }

  while(*ptr != '\0' && isspace(*ptr))
    ptr++;
  
  if (*ptr == '\0' || *ptr == '#') 
    return(0);

  while((token = strtok(ptr," \n,;\t")) != NULL) {
    if (strcmp(token,"-end") == 0)
      return(1);
    ptr = NULL;
    parse_token(token, mode);
  }
  return(0);
}


BOOL
UNCHECKED_PARAM_TABLE::get_bool(const char* name)
{
  return(find(name) != NULL);
}

INT
UNCHECKED_PARAM_TABLE::get_size(const char* name)
{
  PARAM *pe;
  pe = find(name);
  if (pe == NULL) 
    return(0);
  return(pe->size);
}

INT
UNCHECKED_PARAM_TABLE::get_int(const char* name, INT default_value, INT index)
{
  const char* string;
  const char* cp;

  string = get_string(name, NULL, index);
  if (string == NULL) {
    if (default_log)
      fprintf(default_log, "default %s: %d\n", name, default_value);
    return(default_value);
  }
  for(cp = string; *cp != '\0'; cp++) {
    if (isspace(*cp) || *cp == '-' || *cp == '+')
      continue;
    if (isdigit(*cp))
      return(atoi(string));
    else 
      break;
  }
  panic("error: invalid integer for parameter %s: %s\n", name, string);
  return(-1);
}

INT
UNCHECKED_PARAM_TABLE::get_enum(const char* name, INT default_value, const char* *name_table, 
				INT index)
{
  const char* string;
  const char* cp;
  const char* *np = name_table;
  INT i, n;

  string = get_string(name, NULL, index);
  if (string == NULL) {
    for(n=0; np[n] != NULL; n++);
    assert(default_value >= 0 && default_value < n);
    if (default_log)
      fprintf(default_log,
	      "default %s: %s (%d)\n", name, name_table[default_value], default_value);
    return(default_value);
  }
  n = 0;
  while ((cp = *np++) != NULL) {
    if (strcmp(cp, string) == 0)
      return(n);
    n++;
  }

  // try for a number
  cp = string;
  while(isspace(*cp) || *cp == '-' || *cp == '+')
    cp++;
  if (isdigit(*cp)) {
    i = atoi(string);
    if (i >= 0 && i < n) {
      if (default_log)
	fprintf(default_log, 
		"warning: parameter `%s' with number %d instead of name `%s'\n",
		name, i, name_table[i]);
      return(i);
    }
  }

  fprintf(stderr, 
    "PARAM_TABLE ERROR: can not find parameter value `%s' in list:\n", 
	  string);
  for(i=0; name_table[i] != NULL; i++)
    fprintf(stderr, "%s ", name_table[i]);
  fprintf(stderr, "\n");
  panic("PARAM_TABLE::get_enum parameter error");
  return(-1);
}

int*
UNCHECKED_PARAM_TABLE::get_int_array(const char* name)
{
  const char* string;
  int i,n;
  int* array;

  if (!get_bool(name))
    return(NULL);

  n = get_size(name);
  array = new INT[n];
  for(i=0; i<n; i++) {
    string = get_string(name, NULL, i);
    assert(string != NULL);
    array[i] = atoi(string);
  }
  return(array);
}


float *
UNCHECKED_PARAM_TABLE::get_float_array(const char* name)
{
  const char* string;
  int i,n;
  FLOAT* array;

  if (!get_bool(name))
    return(NULL);

  n = get_size(name);
  array = new FLOAT[n];
  for(i=0; i<n; i++) {
    string = get_string(name, NULL, i);
    assert(string != NULL);
    array[i] = atof(string);
  }
  return(array);
}

INT *
UNCHECKED_PARAM_TABLE::get_range_array(const char* name, INT size,
					   INT included, INT excluded)
{
  const char* string;
  int* array;
  int i,j,n;
  int num, prev_num;

  if (!get_bool(name))
    return(NULL);

  array = new INT[size];
  for(i=0; i<size; i++)
    array[i] = excluded;

  n = get_size(name);
  num = -1;
  for(i=0; i<n; i++) {
    string = get_string(name, NULL, i);
    assert(string != NULL);
    prev_num = num;
    num = atoi(string);
    if (num >= 0 && num < size) {
      array[num] = included;
    } else if (num < 0 && num >= -size && prev_num >= 0 && prev_num < -num) {
      num = -num;
      for(j=prev_num; j<num; j++)
	array[j] = included;
    } else
      fprintf(stderr, 
	      "PARAM_TABLE ERROR: bad range value for parameter %s: ... %d %d\n", 
	      name, prev_num, num);
  }
  return(array);
}


INT *
UNCHECKED_PARAM_TABLE::get_mapping(const char* name, INT *map_size,
  INT not_mapped_value)
{
  INT *array;
  INT i,n;
  INT index, value, size;

  if (!get_bool(name))
    return(NULL);

  n = get_size(name);
  size = -1;
  for(i=0; i<n; i += 2) {
    if (get_int(name,-1,i) > size)
      size = get_int(name,-1,i);
  }
  assert(size > 0);
  size++;  // for index 0

  array = new INT[size];
  *map_size = size;
  for(i=0; i<size; i++)
    array[i] = not_mapped_value;

  for(i=0; i<n; i += 2) {
    index = get_int(name,-1,i);
    value = get_int(name,-1,i+1);
    assert(index >= 0 && index < size);
    array[index] = value;
  }

  return(array);
}


float
UNCHECKED_PARAM_TABLE::get_float(const char* name, float default_value, INT index)
{
  const char* string;

  string = get_string(name, NULL, index);
  if (string == NULL) {
    if (default_log)
      fprintf(default_log, "default %s: %f\n", name, default_value);
    return(default_value);
  }
  return(atof(string));
}

const char* 
UNCHECKED_PARAM_TABLE::get_string(const char* name, const char* default_value, INT index)
{
  PARAM *pe;
  pe = find(name);
  if (pe == NULL || index >= pe->size) {
    if (default_value != NULL && default_log)
      fprintf(default_log, "default %s: %s\n", name, default_value);
    return(default_value);
  }
  return(pe->data[index]);
}

CHAR* 
UNCHECKED_PARAM_TABLE::copy_string(const char* name, const char* default_value, INT index)
{
  const char* str;
  str = get_string(name,default_value,index);
  if (str == NULL)
    return(NULL);
  return(strdup(str));
}

void
UNCHECKED_PARAM_TABLE::print(FILE *file, STRING first_char)
{
  PARAM *pe;

  fprintf(file,"%s\n", first_char);
  for(pe = head; pe != NULL; pe = pe->next) {
    print_entry(file, pe, first_char);
  }
  fprintf(file,"%s\n%s\n", first_char, first_char);
  fflush(file);
}

void
UNCHECKED_PARAM_TABLE::print_entry(FILE* file, PARAM* pe, STRING first_char)
{
  const INT maxline = 72;
  INT col = 0;
  INT i;

  fprintf(file, "%s %s = ", first_char, pe->name);
  col = strlen(first_char) + strlen(pe->name)+2;
  for(i=0; i<pe->size; i++) {
    col += strlen(pe->data[i])+1;
    if (col > maxline) {
      fprintf(file,"\n%s   ", first_char);
      col = strlen(pe->data[i])+5;
    }
    fprintf(file, "%s ", pe->data[i]);
  }    
  fprintf(file, "\n");
}

// PARAM_TABLE functions

PARAM_TABLE::PARAM_TABLE(const char* _file_name, const char** param_table)
  : UNCHECKED_PARAM_TABLE(NULL)
{ 
  check_param = 0;
  if (param_table) {
    master_parameter_list.parse_cmd_line(get_num_parameters(param_table)*2-1, 
				param_table, OVERRIDE);
    check_param = 1;
  }
  if (_file_name) {
    file_name = strdup(_file_name);
    parse_file(file_name);
  }
}

PARAM* 
PARAM_TABLE::put(const char* name) { 
  if (check_param && !master_parameter_list.get_bool(name))  {
    printf("** PARAM_TABLE ERROR: unknown parameter: %s\n", name);
    if (file_name)
      printf("pfile = %s\n", file_name);
    help_message();
    exit(-1);
  }
  return(UNCHECKED_PARAM_TABLE::put(name));
}


PARAM* 
PARAM_TABLE::find(const char* name) 
{ 
  if (check_param && !master_parameter_list.get_bool(name))  {
    printf(
	   "** PARAM_TABLE ERROR: unregistered param used: %s\n", 
	   name);
    if (file_name)
      printf("param_file = %s\n", file_name);
    printf("Use noway -help for a listing of all params\n");
    exit(1);
  }
  return(UNCHECKED_PARAM_TABLE::find(name));
}


void 
PARAM_TABLE::help_message() {
  printf("\n\n");
  printf("****************************\n");
  printf("*** List of parameters:  ***\n");
  printf("****************************\n");
  master_parameter_list.print(stdout);
  printf("****************************\n\n");
}


INT 
PARAM_TABLE::get_num_parameters(const char** p)
{
  INT i;
  for(i=0; p[i] != NULL; i += 2);
  return(i/2);
}


const char* itoa(INT i)
{
  static char ret [30];
  sprintf(ret,"%d",i);
  return(ret);
}

