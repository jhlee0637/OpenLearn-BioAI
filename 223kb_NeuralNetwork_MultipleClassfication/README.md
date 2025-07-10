[한국어](./README_kor.md)
# Neural Network Multiple Classification - Fetal Health

## Overview
This project implements a PyTorch-based neural network for multi-class classification of fetal health using Cardiotocogram (CTG) data. The model classifies fetal health into three categories: Normal, Suspect, and Pathological.

## Dataset
- **Source**: [Fetal Health Classification - Kaggle](https://www.kaggle.com/datasets/andrewmvd/fetal-health-classification)
- **Size**: 223KB
- **Description**: Cardiotocogram exams classified by three expert obstetricians into 3 classes
- **Classes**: 
  - Normal
  - Suspect  
  - Pathological

## Model Architecture
- **Model Name**: FetalNet
- **Type**: Multi-layer Neural Network (PyTorch)
- **Task**: Multi-class Classification
- **Features**: Standardized using StandardScaler
- **Train/Validation Split**: 80/20

## Environment & Dependencies

### System Specifications
- **OS**: Rocky Linux
- **CPU**: Intel Core i5 11th Gen 1135G7 (2.40GHz)
- **Memory**: 8 GB
- **Graphics**: Intel Iris Xe Graphics

### Python Environment
| Package | Version |
|:-------:|:-------:|
| Python | 3.12.10 |
| PyTorch | 2.7.1+cpu |
| NumPy | 2.3.0 |
| Pandas | 2.3.0 |
| Scikit-learn | 1.7.0 |
| KaggleHub | 0.3.12 |

## Files
- `pytorch_model.ipynb`: Main Jupyter notebook containing the complete implementation
- `README.md`: This documentation file

## Usage

### 1. Run the Notebook
```bash
jupyter notebook pytorch_model.ipynb
```

### 2. Dataset Download
The notebook automatically downloads the dataset using KaggleHub:
```python
path = kagglehub.dataset_download("andrewmvd/fetal-health-classification")
```

## Implementation Steps
1. **Data Loading**: Download and load fetal health dataset
2. **Data Preprocessing**: 
   - Feature scaling using StandardScaler
   - Train/validation split
3. **Model Definition**: FetalNet neural network architecture
4. **Training**: Model training with optimization
5. **Evaluation**: Performance assessment on validation set

## Key Features
- **Data Preprocessing**: Automated feature scaling and data splitting
- **PyTorch Implementation**: Modern deep learning framework
- **Multi-class Classification**: Handles 3-class fetal health prediction
- **Reproducible Results**: Structured notebook with clear documentation

## Practice Reference
**Date**: June 18, 2025  
**Instructor**: ChoiTaeOn
