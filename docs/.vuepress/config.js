// .vuepress/config.js
module.exports = {
	locales: {
		'/': {
			lang: 'en-US',
			title: 'ReacNetGenerator',
			description: 'An automatic generator of reaction network for reactive molecular dynamics simulation'
		},
		'/zh/': {
			lang: 'zh-CN',
			title: 'ReacNetGenerator',
			description: '反应动力学模拟的反应网络自动生成器'
		}
	},
	head: [
    ['link', { rel: 'icon', href: `/reacnetgen.svg` }],
    ['meta', { name: 'apple-mobile-web-app-capable', content: 'yes' }],
    ['link', { rel: 'apple-touch-icon', href: `/reacnetgen.svg` }],
    ['meta', { name: 'msapplication-TileImage', content: '/reacnetgen.svg' }]
  ],
	themeConfig: {
		locales: {
			'/': {
				selectText: 'Languages',
				label: 'English',
				nav: [
					{ text: 'Home', link: '/' },
					{ text: 'Report', link: '/r.html' },
					{ text: 'Article', link: 'https://doi.org/10.26434/chemrxiv.7421534'},
					{ text: 'Group', link: 'http://computechem.cn'},
				]
			},
			'/zh/': {
				selectText: '语言',
				label: '中文',
				nav: [
					{ text: '主页', link: '/zh/' },
					{ text: '分析结果', link: '/r.html' },
					{ text: '论文', link: 'https://doi.org/10.26434/chemrxiv.7421534'},
					{ text: '课题组', link: 'http://computechem.cn'},
				]
			},
		}
	},
}
