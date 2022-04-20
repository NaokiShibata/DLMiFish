# DLMiFish
NCBIから真骨類 (Teleostei, taxid=[32443](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=32443&lvl=3&lin=f&keep=1&srchmode=1&unlock))、軟骨魚類 (Chondrichthyes, taxid=[7777](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=7777&lvl=3&lin=f&keep=1&srchmode=1&unlock))、円口類(Cyclostomata, taxid=[1476529](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=1476529&lvl=3&lin=f&keep=1&srchmode=1&unlock))の12s rRNA領域の配列データをダウンロードするパッケージです。`Biopython`に搭載されている`Entrez`の各関数を使って、上記3分類群の配列を取得した後、Genbankの情報をもとに自動で12s rRNA領域を切り出し抽出しています。また、Option機能では、`Cutadapt`によるプライマー配列の一部削除(データ軽量化のため)と`seqkit`による配列データの配列長を視覚化したりもできます。

windows10とUbuntu 20.04 LTSでの動作を確認しています。

## 動作環境や環境整備に必要なもの(抜けが確認出来れば随時更新します)

- python3.6 以上 (f string使用のため)
  - mamba : 最新版 or pythonのバージョンにあったもの
  - pandas : 最新版 or pythonのバージョンにあったもの
  - biopython : 最新版 or pythonのバージョンにあったもの
  - tqdm : 最新版 or pythonのバージョンにあったもの
  - 
- seqkit version2.0以上
  ### Install
  1. Homeデレクトリに移動
  ```bash
  cd
  ```
  2. 最新版をダウンロード(要version確認)
  ```bash
  # 2022/4/21
  # Windows 10
  wget https://github.com/shenwei356/seqkit/releases/download/v2.2.0/seqkit_linux_arm64.tar.gz
  ```
  3. 解凍
  ```
  tar -zxvf *arm64.tar.gz
  ```
  4. プログラムの置き場に移動
  ```
  sudo mv seqkit ~/usr/local/bin/
  ```
  `conda`が使える方は以下で一発
  ```bash
  conda install -c bioconda seqkit -y
  ```
  `mamba`でもよい(むしろパッケージ管理にはこちらがおすすめ
  mamba install -c bioconda seqkit -y
  ```
  5. 確認
  ```
  seqkit version
  # seqkit v2.0.0
  ```
  
  
  
  
