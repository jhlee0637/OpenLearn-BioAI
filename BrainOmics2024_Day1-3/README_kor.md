## BrainOmics2024: Single-cell & Spatial Transcriptomics 분석 실습
### Abstract
이 레포지토리는 bioinformatics 스터디 그룹의 활동 일환으로 작성되었습니다.    
2024년 이탈리아 밀라노의 Human Technopole에서 진행된 [Brain Omics 2024 Course](https://github.com/BrainOmicsCourse/BrainOmics2024)의 실습 내용을 바탕으로, 해당 코스의 내용을 재구성하고 구현하는 것에 초첨을 맞추고 있습니다.    

### 목표
- 코드 주석을 통한 scRNA-seq 분석 파이프라인 완전 이해
- 공간 전사체학 데이터 분석 방법론 학습
- pandas, numpy, scanpy를 활용한 생물정보학 데이터 처리
- PCA, UMAP을 이용한 차원 축소 및 고급 시각화 기법

### Day 0: 환경 설정
- Conda 환경 설정
- AnnData 객체 최적화
- Jupyter Notebook 설정

### Day 1: Single-cell RNA-seq 기초 분석
- Count matrix에서 클러스터까지의 전체 워크플로우
- 데이터 필터링, 정규화, 배치 효과 제거
- 클러스터 특성화 및 마커 유전자 발견

### Day 2: 고급 Single-cell 분석
- Multiplet 검출 및 제거
- Metacell 생성
- 궤적 추론 (Trajectory inference)
- 전사인자 활성 예측
