/*
 * Mole - Mobile Organic Localisation Engine
 * Copyright 2010 Nokia Corporation.
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 * http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "completer.h"

MultiCompleter::MultiCompleter
(QObject *parent)
     : QCompleter(parent)
 {

   //sep = QLatin1String(",");
   sep = QLatin1String("");
 }

 MultiCompleter::MultiCompleter(QAbstractItemModel *model, QObject *parent)
     : QCompleter(model, parent)
 {
 }

void MultiCompleter::setSeparator (QLatin1String _sep) {
  sep = _sep;
}

QString MultiCompleter::separator() const
 {
     return sep;
 }

 QStringList MultiCompleter::splitPath(const QString &path) const
 {
     if (sep.isNull()) {
         return QCompleter::splitPath(path);
     }

     return path.split(sep);
 }

QString MultiCompleter::pathFromIndex(const QModelIndex &index) const
 {
     if (sep.isNull()) {
         return QCompleter::pathFromIndex(index);
     }

     // navigate up and accumulate data
     QStringList dataList;
     for (QModelIndex i = index; i.isValid(); i = i.parent()) {
         dataList.prepend(model()->data(i, completionRole()).toString());
     }

     return dataList.join(sep);
 }


 void MultiCompleter::setSeparator(const QString &separator)
 {
     sep = separator;
 }

