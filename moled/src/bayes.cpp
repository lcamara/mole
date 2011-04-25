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

#include "bayes.h"


qreal Histogram::findBucket (qreal value) {
  return (floor(value / bucketWidth)) * bucketWidth;
}

void Histogram::add(Histogram *h) {
  Q_ASSERT (bucketWidth == h->bucketWidth);

  for (QMap<qreal,int>::const_iterator iter = h->m_value2count.constBegin(); iter != h->m_value2count.constEnd(); ++iter) {
    //add (iter->first(), iter->second(), false);
    Q_ASSERT(m_size >= 0);
  }
}

void Histogram::add (qreal value, int count, bool doCoerce) {
  // to get more exact mean, remember values before we coerce them
  m_sum += (value * count);

  if (doCoerce) {
    value = findBucket (value);
  }

  QMap<qreal,int>::iterator iter = m_value2count.find(value);

  if (iter == m_value2count.end()) {
    //m_value2count.insert (QPair<qreal,int>(value,count));
    m_values.push_back (value);
  } else {
    //iter->second += count;
  }
  m_size += count;

}

qreal Histogram::exactMean () const
{
  if (m_size == 0) return 0.;
  return (m_sum / m_size);
}

qreal Histogram::mean() const
{
  qreal sum = 0.0;
  qreal value;
  int count;
  for (QMap<qreal, int>::const_iterator it = m_value2count.constBegin(); it != m_value2count.constEnd(); ++it) {
    //value = it->first;
    //count = it->second;
    sum += value * count;
  }
  return sum / m_size;
}

qreal Histogram::var() const
{
  if (m_size == 1) {
    return 0;
  }

  qreal meanValue = mean();
  qreal squareSum = 0.0;
  qreal value;
  qreal deviation;
  int count;

  for (QMap<qreal, int>::const_iterator it = m_value2count.constBegin(); it != m_value2count.constEnd(); ++it) {
    //value = it->first;
    //count = it->second;
    deviation = value - meanValue;
    squareSum += deviation * deviation * count;
  }

  qreal varValue = squareSum / (m_size-1);

  Q_ASSERT(varValue >= 0); // Otherwise, something is wrong...

  return varValue;
}

qreal Histogram::std() const
{
  return sqrt(var());
}

qreal Histogram::percentile(qreal pct) {
  qSort (m_values.begin(), m_values.end());
  int currentIndex = 0;
  int targetIndex = int(m_size * pct);
  double currentValue;

  if (m_size == 0)
      return 0;

  Q_ASSERT (m_values.size() > 0);

  for (QMap<qreal,int>::const_iterator iter = m_value2count.constBegin(); iter != m_value2count.constEnd(); ++iter) {
    //currentIndex += iter->second;
    //currentValue = iter->first;

    if (currentIndex >= targetIndex)
      break;
  }
  return currentValue;
}

// Frequency of a value.
qreal Histogram::getFrequency(qreal value)
{
  QMap<qreal, int>::const_iterator it = m_value2count.constFind(value);

  if (it != m_value2count.constEnd()) {
    int ct = it.value();
    return ct * 1.0;
  } else {
    return 0.0;
  }
}

// Probability of a value. This is just a maximum-likelihood estimate
// (= count/total_size). If want to have a kernel histogram, do something
// else.
qreal Histogram::probability(qreal value)
{
  QMap<qreal, int>::const_iterator it = m_value2count.constFind(value);

  if (it != m_value2count.constEnd()) {
    int ct = it.value();
    return (qreal)ct/m_size;
  } else {
    return 0.0;
  }
}

/********** Kernel Signature **********/
/*void KernelHistogram::add(int value)
{
  Histogram::add(value);
  m_kernelUpdated = false;
}*/

qreal KernelHistogram::getFrequency(int value)
{
  int histSize = Histogram::getSize();
  qreal freq = 0.0;

  if (m_kernelUpdated) {
    freq = m_kernelPdf[value] / m_kernelPdfWeight * histSize;

  } else {			// Update kernel estimate first.
    m_kernelPdf.clear();
    m_kernelPdfWeight = 0.0;

    QVector<qreal> histKeys = getKeyValues();

    for (QVector<qreal>::iterator it = histKeys.begin(); it != histKeys.end(); ++it) {
      qreal key = *it;
      int keyInteger = floor(key+.5); // Exact integer key by rounding... hopefully.
      for (uint dev = 0; dev < m_kernelHalfWidth; ++dev) {
        if (dev != 0) {
          // It's assumed that the default constructor of "qreal" will fill in 0.0
          int devKey_p = keyInteger + dev;
          qreal kerValue = m_defaultKernel[dev] * Histogram::getFrequency(key);
          m_kernelPdf[devKey_p] += kerValue;
          m_kernelPdfWeight += kerValue;

          int devKey_m = keyInteger - dev;
          kerValue = m_defaultKernel[dev] * Histogram::getFrequency(key);
          m_kernelPdf[devKey_m] += kerValue;
          m_kernelPdfWeight += kerValue;
        }
        else {
          qreal kerValue = m_defaultKernel[0] * Histogram::getFrequency(key);
          m_kernelPdf[keyInteger] += kerValue;
          m_kernelPdfWeight += kerValue;
        }
      }
    }
    m_kernelUpdated = true;

    freq = m_kernelPdf[value] / m_kernelPdfWeight * histSize;
  }

  return freq;
}

qreal KernelHistogram::m_defaultKernel[6] =
              {0.2006, 0.1770, 0.1217, 0.06511, 0.02714, 0.008812};
//                0        +-1     +-2      +-3     +-4      +-5
uint KernelHistogram::m_kernelHalfWidth = 6;

BayesSignature::BayesSignature(const AP_Scan *scan)
  : m_numberOfImplicatedScans(0)
  , m_numberOfTotalReadings(0)
{
  addScan(scan->readings_map);
}

BayesSignature::BayesSignature()
  : m_numberOfImplicatedScans(0)
  , m_numberOfTotalReadings(0)
{
  qDebug() << "BayesSignature constructor";
  //addScan(scan->readings_map);
}

BayesSignature::~BayesSignature()
{
  clear();
}

void BayesSignature::clear()
{
  m_numberOfImplicatedScans = 0;
  m_numberOfTotalReadings = 0;
  m_macs.clear();

  // Delete histograms.
  qDeleteAll(m_mac2hist.begin(), m_mac2hist.end());
  m_mac2hist.clear();
}

Histogram* BayesSignature::histogram(QString apMac)
{
  QMap<QString,Histogram*>::const_iterator it = m_mac2hist.find(apMac);

  if(it != m_mac2hist.end())
    return it.value();

  return 0;
}

void BayesSignature::addScan(QMap<QString,quint32> *readingMaps)
{
  QMapIterator<QString,quint32> it (*readingMaps);
  QMap<QString,Histogram*>::const_iterator itHist;

  while(it.hasNext()) {
    it.next();
    QString apMac = it.key();
    quint32 rssi = it.value();
    m_macs.insert(apMac);
    itHist = m_mac2hist.find(apMac);

    //Find or create a histogram
    Histogram *pHist;
    if(itHist == m_mac2hist.end()) {  //apMac is not found
      pHist = newHistogram();
      m_mac2hist.insert(apMac, pHist);
    } else {
      pHist = itHist.value();
    }
    //Update the histogram
    pHist->add(rssi);

    //Recording the total number of readings.
    ++m_numberOfTotalReadings;
  }

  //Recording the total number of scans.
  ++m_numberOfImplicatedScans;
}
