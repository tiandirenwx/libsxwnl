import { hapTasks } from '@ohos/hvigor-ohos-plugin';
import { execSync } from 'child_process';
import * as path from 'path';

// 构建前链接仓库根 assets/bazi/ 共享资源 (与 Android/iOS 同源)
const repoRoot = path.resolve(__dirname, '../..');
execSync(`bash "${path.join(repoRoot, 'scripts/sync_bazi_assets.sh')}" harmony`, {
  stdio: 'inherit',
});

export default {
  system: hapTasks,
};
