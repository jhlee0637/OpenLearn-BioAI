[English](./README.md)
---

# 신경망 다중 분류 - 태아 건강 예측

## 개요
이 프로젝트는 심전도(CTG) 데이터를 사용하여 태아 건강을 다중 분류하는 PyTorch 기반 신경망을 구현합니다. 모델은 태아 건강을 정상, 의심, 병리학적 세 가지 범주로 분류합니다.

## 데이터셋
- **출처**: [Fetal Health Classification - Kaggle](https://www.kaggle.com/datasets/andrewmvd/fetal-health-classification)
- **크기**: 223KB
- **설명**: 3명의 전문 산부인과 의사가 3개 클래스로 분류한 심전도 검사 결과
- **분류**: 
  - 정상 (Normal)
  - 의심 (Suspect)
  - 병리학적 (Pathological)

## 모델 구조
- **모델명**: FetalNet
- **유형**: 다층 신경망 (PyTorch)
- **작업**: 다중 클래스 분류
- **특성**: StandardScaler를 사용한 표준화
- **훈련/검증 분할**: 80/20

## 환경 및 의존성

### 시스템 사양
- **운영체제**: Rocky Linux
- **CPU**: Intel Core i5 11th Gen 1135G7 (2.40GHz)
- **메모리**: 8 GB
- **그래픽**: Intel Iris Xe Graphics

### Python 환경
| 패키지 | 버전 |
|:-------:|:-------:|
| Python | 3.12.10 |
| PyTorch | 2.7.1+cpu |
| NumPy | 2.3.0 |
| Pandas | 2.3.0 |
| Scikit-learn | 1.7.0 |
| KaggleHub | 0.3.12 |

## 파일 구성
- `pytorch_model.ipynb`: 전체 구현이 포함된 메인 Jupyter 노트북
- `README.md`: 영문 문서
- `README_kor.md`: 한글 문서 (현재 파일)

## 사용법

### 1. 노트북 실행
```bash
jupyter notebook pytorch_model.ipynb
```

### 2. 데이터셋 다운로드
노트북에서 KaggleHub를 사용하여 자동으로 데이터셋을 다운로드합니다:
```python
path = kagglehub.dataset_download("andrewmvd/fetal-health-classification")
```

## 구현 단계
1. **데이터 로딩**: 태아 건강 데이터셋 다운로드 및 로드
2. **데이터 전처리**: 
   - StandardScaler를 사용한 특성 스케일링
   - 훈련/검증 데이터 분할
3. **모델 정의**: FetalNet 신경망 아키텍처 구성
4. **훈련**: 최적화를 통한 모델 훈련
5. **평가**: 검증 세트에서 성능 평가

## 주요 특징
- **데이터 전처리**: 자동화된 특성 스케일링 및 데이터 분할
- **PyTorch 구현**: 최신 딥러닝 프레임워크 사용
- **다중 클래스 분류**: 3개 클래스 태아 건강 예측 처리
- **재현 가능한 결과**: 명확한 문서화가 포함된 구조화된 노트북

## 기술적 세부사항

### 데이터 전처리
- **특성 스케일링**: 모든 특성을 표준화하여 모델 성능 향상
- **데이터 분할**: 80%는 훈련용, 20%는 검증용으로 분할
- **텐서 변환**: PyTorch 텐서로 데이터 타입 변환

### 모델 아키텍처
- **입력층**: 특성 수에 맞는 입력 노드
- **은닉층**: 다층 완전 연결 신경망
- **출력층**: 3개 클래스에 대한 소프트맥스 출력
- **손실 함수**: 교차 엔트로피 손실
- **최적화**: Adam 옵티마이저 사용

### 성능 지표
- **정확도**: 전체 예측 정확도
- **클래스별 성능**: 각 클래스에 대한 정밀도, 재현율, F1-점수
- **혼동 행렬**: 분류 성능의 시각적 표현

## 학습 목표
이 프로젝트를 통해 다음을 학습할 수 있습니다:
- PyTorch를 사용한 신경망 구현
- 의료 데이터를 활용한 다중 분류 문제 해결
- 데이터 전처리 및 특성 엔지니어링
- 모델 훈련 및 평가 과정
- Jupyter 노트북을 활용한 데이터 과학 워크플로우

## 향후 개선 방향
- **하이퍼파라미터 튜닝**: 그리드 서치 또는 베이지안 최적화
- **교차 검증**: K-fold 교차 검증을 통한 더 robust한 평가
- **앙상블 방법**: 여러 모델을 결합하여 성능 향상
- **시각화 개선**: 더 자세한 결과 분석 및 시각화
- **배포**: 웹 애플리케이션 또는 API로 모델 배포

## 참고 자료
- [PyTorch 공식 문서](https://pytorch.org/docs/)
- [Scikit-learn 사용자 가이드](https://scikit-learn.org/stable/user_guide.html)
- [Kaggle 데이터셋](https://www.kaggle.com/datasets/andrewmvd/fetal-health-classification)

## 실습 정보
**날짜**: 2025년 6월 18일  
**강사**: ChoiTaeOn
