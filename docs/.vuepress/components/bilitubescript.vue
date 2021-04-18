<template>
  <component v-if="VueScriptComponent" :is="VueScriptComponent" script="<script src='https://cdn.jsdelivr.net/npm/bilitube@0/dist/bilitube.min.js' defer></script>"/>
</template>
<script>
export default {
  data() {
    return {
      VueScriptComponent: null
    }
  },
  mounted () {
    import('vue-script-component').then(module => {
      this.VueScriptComponent = module.default
    })
  }
}
</script>
