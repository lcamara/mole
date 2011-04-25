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

#ifndef BAYES_H_
#define BAYES_H_

#include "scan.h"

#include <QtCore>

class Histogram
{

public:
  const qreal bucketWidth;
  Histogram (qreal bw) : bucketWidth(bw) { m_size = 0; m_sum = 0.0;}

  void add(quint32 value) { add((qreal)value,1); }
  void clear() { m_value2count.clear(); m_size = 0; m_values.clear(); }
  qreal exactMean() const;      // Mean before coerced.
  qreal mean() const;           // Mean after binned.
  qreal var() const;            // Sample variance.
  qreal std() const;            // Sample standard deviation.
  qreal percentile(qreal pct);
  int getSize() const { return m_size; }
  qreal getFrequency(qreal value);
  virtual qreal getFrequency(int value) const { return getFrequency((qreal)value); }
  qreal probability(qreal value);
  qreal probability(int value) const { return probability((qreal)value); }
  QVector<qreal> getKeyValues() const { return m_values; }

private:
  int m_size;
  qreal m_sum;

  // TODO some kind of sorted list would be better
  QVector<qreal> m_values;
  QMap<qreal,int> m_value2count;
  void add(qreal value, int count, bool doCoerce);
  qreal findBucket(qreal value);

  void add(Histogram *);
  void add(qreal value, int count) { add(value, count, true); }
  void add(qreal value) { add(value, 1); }
  void add(int value, int count) { add((qreal)value,count); }
};


// Note that an output of get_freq() is multiplied by the total number
// of readings so that the total weight of this kernel estimate
// function matches to the original Histogram.
//
// Kernel estimation is done is lazy manner.
class KernelHistogram : public Histogram
{

public:
  KernelHistogram(qreal bw) : Histogram(bw), m_kernelPdfWeight(0.0), m_kernelUpdated(false) { }
  ~KernelHistogram() { }

  //void add(int value);              // Add a value and tag kernel_updated as false.
  qreal getFrequency(int value);   // Returns prob * weight;

private:
  QMap<int, qreal> m_kernelPdf;
  qreal m_kernelPdfWeight;
  bool m_kernelUpdated;

  // Kernel.
  static qreal m_defaultKernel[6];
  static uint m_kernelHalfWidth;
};

class BayesSignature
{

public:
  BayesSignature(const AP_Scan *scan);
  BayesSignature();
  ~BayesSignature();

  // Add a scan to this signature.
  void addScan(QMap<QString,quint32> *readingMaps);

  void clear();

  // Return the histogram for mac.
  Histogram* histogram(QString apMac);

  Histogram* newHistogram() { return new Histogram(1.); // XXX Need floating point histogram?
                            }
  QSet<QString> macs() const { return m_macs; }
  QMap<QString,Histogram*> mac2hist() const { return m_mac2hist; }

private:
  QSet<QString> m_macs;
  QMap<QString,Histogram*> m_mac2hist;
  int m_numberOfImplicatedScans;
  quint32 m_numberOfTotalReadings;
};

#endif /* BAYES_H_ */
