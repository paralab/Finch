/*
  Copyright 2014-2016 Baskar Ganapathysubramanian

  This file is part of TALYFem.

  TALYFem is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as
  published by the Free Software Foundation, either version 2.1 of the
  License, or (at your option) any later version.

  TALYFem is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with TALYFem.  If not, see <http://www.gnu.org/licenses/>.
*/
// --- end license text --- //
#include <talyfem/common/pack_comm.h>

namespace TALYFEMLIB {

void getLine(FILE* fp, char* str) {
  char tmp;
  int i;
  for (i = 0;; i++) {
    if (!feof(fp)) {
      int nread = static_cast<int>(fread(&tmp, sizeof(char), 1, fp));
      if (nread != 1)
        break;
    } else {
      break;
    }
    if (tmp != '\n') {
      str[i] = tmp;
    } else {
      break;
    }
  }
  str[i] = '\0';
}

char* getParameter(FILE* fp, const char* varname, char* var,
                   const char* divideChar) {
  var[0] = '\0';
  char line[1024];
  for (;;) {
    getLine(fp, line);
    char* pHead = strstr(line, varname);
    if (pHead) {
      pHead = strstr(line + strlen(varname), divideChar) + 1;
      trimSpace(pHead, var);
      break;
    }
  }
  return var;
}

void trimSpace(char* str, char* s) {
  int end = static_cast<int>(strlen(str));
  int start = 0;
  int i, j;
  // trim right
  for (i = end - 1; i >= 0; i--) {
    if (str[i] == ' ') {
      end--;
    } else {
      break;
    }
  }
  // trim left
  for (i = 0; i < end; i++) {
    if (str[i] == ' ') {
      start++;
    } else {
      break;
    }
  }
  // copy
  for (i = start, j = 0; i < end; i++, j++) {
    s[j] = str[i];
  }
  s[j] = '\0';
}

int divideParameters(char* str, char* position[], const char* divideChar) {
  int len = static_cast<int>(strlen(str));
  int i;
  int n = 0;
  int wait = 1;
  for (i = 0; i < len; i++) {
    const char* p = strchr(divideChar, str[i]);
    if (p) {
      wait = 1;
      str[i] = '\0';
    } else if (wait) {
      position[n] = str + i;
      n++;
      wait = 0;
    }
  }
  // ~ *position[n]='\0'; //position[n]='\0'; ->
  // warning: expression which evaluates to zero treated as a null
  // pointer constant of type 'char *'
  position[n] = NULL;
  return n;
}

}  // namespace TALYFEMLIB
