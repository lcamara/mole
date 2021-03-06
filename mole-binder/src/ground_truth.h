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

#ifndef GROUND_TRUTH_H
#define GROUND_TRUTH_H

#include <QtGui>
#include <QtDBus>
#include <QWidget>

#include "common.h"
#include "settings.h"

class GroundTruth : public QWidget
{
    Q_OBJECT
public:
  explicit GroundTruth(QWidget *parent = 0, QString groundTruthFilename = 0);
    ~GroundTruth ();

signals:

public slots:
  void handle_quit();

private:

  QList<QPushButton *> buttons;
  void build_ui (QFile &file);

private slots:
  void buttonClicked();

};

#endif // GROUND_TRUTH_H
